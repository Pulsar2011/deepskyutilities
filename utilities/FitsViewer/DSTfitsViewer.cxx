#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <map>
#include <set>

#ifdef __WITH_QT__

#include <QtConcurrent/QtConcurrent>
#include <DSTdisplay/QtFITSviewer.h>
#include <QMetaObject>
#include <QMainWindow>
#include <QStatusBar>
#include <QSlider>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QProxyStyle>
#include <QPainter>
#include <QStyleOption>
#include <QPolygon>
#include <QSplitter>
#include <QComboBox>
#include <QLineEdit>
#include <QListWidget>
#include <QClipboard>
#include <QMenu>
#include <QGuiApplication>
#include <QToolButton>
#include <QKeyEvent>
#include <QRadioButton>
#include <QAbstractTableModel>
#include <QTableView>
#include <QFutureWatcher>
#include <QIcon>
#include <QPixmap>
#include <QPushButton>
#include <QCheckBox>
#include <QTableWidget>
#include <QHeaderView>
#include <QFormLayout>
#include <QToolBar>
#include <QTabWidget>
#include <QMessageBox>
#include <QSpinBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QGroupBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QDockWidget>
#include <QGridLayout>
#include <QScrollArea>
#include <QStackedWidget>
#include <QSignalBlocker>
#include <QStringList>
#include <QFileInfo>
#include <QDoubleValidator>
#include <QLocale>
#include <QFontDatabase>
#include <QPointer>
#include <DSTdisplay/QCustomContours.h>
#include <DSTdisplay/QtFitsPloting.h>

#endif

#include <limits>

#include <DSTCore/OptionManager.h>
#include <DSTCore/DSS_version.h>
#include <DSTools/DST_version.h>
#include <DSTdisplay/DSD_version.h>
#include <DSTfits/DSF_version.h>
#include <DSTfits/FITSmanager.h>
#include <DSTfits/FITSimg.h>
#include <DSTfits/FITShdu.h>
#include <DSTfits/FITSnumeric.h>
#include <DSTmath/DSTmath.h>

#ifdef __WITH_ASDF__
#include <cstring>
#include <memory>
#include <asdf2fits/ASDFmanager.h>
#endif


#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future> // added for std::async

typedef std::thread THREAD;

std::thread _draw_histo;
std::thread _draw_grid;
std::thread _print_info;
std::thread _print_;

// Set to true before joining background threads when switching images so that
// DrawWCSgrid can abort its expensive work rather than running to completion.
static std::atomic<bool> _cancel_wcs{false};

using namespace DSL;

#ifdef __WITH_ASDF__

// Opens a FITS or an ASDF file and hands out the same DeepSkyFits objects for both. The format is
// taken from the file content, not its name; anything that isn't ASDF (including the CFITSIO
// extended file names, which can't be opened as a plain path) goes to FITSmanager.
class DataFile
{
    std::unique_ptr<FITSmanager> fits;
    std::unique_ptr<ASDFmanager> asdf;

    static bool isASDF(const std::string& path)
    {
        char magic[5] = {};
        std::ifstream in(path, std::ios::binary);
        in.read(magic, sizeof(magic));
        return in.gcount() == sizeof(magic) && std::memcmp(magic, "#ASDF", sizeof(magic)) == 0;
    }

public:
    explicit DataFile(const std::string& path)
    {
        if (isASDF(path))
            asdf = std::make_unique<ASDFmanager>(path);
        else
            fits = std::make_unique<FITSmanager>(path);
    }

    FITSmanager* fitsFile() { return fits.get(); }      //!< nullptr when the file is ASDF

    int NumberOfHeader() { return asdf ? asdf->NumberOfHeader() : fits->NumberOfHeader(); }

    std::shared_ptr<FITShdu> GetHeaderAtIndex(const int& i)
    {
        if (asdf) return asdf->GetHeaderAtIndex(i);
        return fits->GetHeaderAtIndex(i);
    }

    std::shared_ptr<FITScube> GetImageAtIndex(const int& i)
    {
        if (asdf) return asdf->GetImageAtIndex(i);
        return fits->GetImageAtIndex(i);
    }

    std::shared_ptr<FITStable> GetTableAtIndex(const int& i)
    {
        if (asdf) return asdf->GetTableAtIndex(i);
        return fits->GetTableAtIndex(i);
    }

    void Close()
    {
        if (asdf) asdf->Close();
        else      fits->Close();
    }
};

#else

using DataFile = FITSmanager;

#endif

typedef const std::vector<std::string> flist;
typedef std::map<char,std::pair<size_t,size_t> > roi;

size_t histo_bins=100;
double minZscale=std::numeric_limits<double>::quiet_NaN();
double maxZscale=std::numeric_limits<double>::quiet_NaN();

const std::string DSS_need = "v2.0.0";
const std::string DSD_need = "v1.2.4";
const std::string DSF_need = "v1.5.0";
const std::string ThisVersion = "v2.2.1";

class UOption : public DOption
{
public:
    
    enum rmode
    {
        NONE = 0x0,
        SUM  = 0x1,
        AVG  = 0x3,
    };
    
    enum drawing_mode
    {
        LINEAR=0x0,
        OFF   =0xF,
    };
    
private:
    int16_t HDU;
    size_t RBIN;
    rmode MODE;
    bool do_Print;
    bool a_info;
    bool hdu_explicit;

public:
    
    static drawing_mode ddraw;
    roi ROI;
    
    UOption():DOption(),HDU(1),RBIN(0),MODE(rmode::NONE), do_Print(false),a_info(false),hdu_explicit(false),ROI()
    {
#ifdef __WITH_ASDF__
        AddDescription("Open a FITS or ASDF header block and display data.");
#else
        AddDescription("Open a FITS header block and display data.");
#endif
        AddOption("--hdu",DParam::fparameter::uint8,"Header block to be open. If not given, the primary header will be opened");
        AddOption("--bins/-b",DParam::fparameter::uint16,"Number of bin to display the histograms of the images.");
        AddOption("--roi",{DParam::fparameter::xoffset,DParam::fparameter::xsize,DParam::fparameter::yoffset,DParam::fparameter::ysize},"Define a refgion of interest to be ploted. Other data outside this ROI are disregarded. The ROI shall be given in pixel not in WCS");
        AddOption("--print","Print pixel content to stdout.");
        AddOption("--info","Print statistical information on the array distribution.");
        AddOption("--draw","Disable drawing!");
        AddOption("--zscale/-z",{DParam::fparameter::min,DParam::fparameter::max},"Select the min and max value for the color scale.");
#ifdef __WITH_ASDF__
        AddParam("Fits_Filename","FITS file URL and name, or ASDF file name, that will be opened for analyses.");
#else
        AddParam("Fits_Filename","FITS file URL and name that will be opened for analyses.");
#endif
        
        DOption::prg_name += "FitsViewer -- "+ThisVersion;

        if(DSS::gGIT::version() < DSS_need)
            throw std::runtime_error("[ERROR] UOption::UOption: Incompatible DSS version, please update your DSS installation to "+DSS_need+" at minimum.");
        if(DSD::gGIT::version() < DSD_need)
            throw std::runtime_error("[ERROR] UOption::UOption: Incompatible DSD version, please update your DSD installation to "+DSS_need+" at minimum.");
        if(DSF::gGIT::version() < DSF_need)
            throw std::runtime_error("[ERROR] UOption::UOption: Incompatible DSF version, please update your DSF installation to "+DSF_need+" at minimum.");
    }
    
    virtual ~UOption(){};
    
    bool UserParameters(std::vector<std::string>::const_iterator& it, const std::vector<std::string>& arg)
    {
        
        if(!it->compare("--hdu"))
        {
            if(std::next(it) != arg.cend()){it++;}else{return false;}
            HDU = std::stoi(*it);
            hdu_explicit = true;
            std::cout<<"   \033[34m• Open header block:\033[0m"<<HDU<<std::endl;
        }
        else if(!it->compare("--info") ||
                !it->compare("-i"))
        {
            a_info |= true;
        }
        else if(!it->compare("--bin") || !it->compare("--bins") || !it->compare("-b"))
        {
            if(std::next(it) != arg.cend()){it++;}else{return false;}
            histo_bins = std::stoi(*it);
        }
        else if(!it->compare("--draw"))
        {
            UOption::ddraw = UOption::drawing_mode::OFF;
        }
        else if(!it->compare("--zscale") || !it->compare("-z"))
        {
            if(std::next(it) != arg.cend()){it++;}else{return false;}
            if(it->find_first_not_of("+-0123456789.eE") != std::string::npos)
                return false;
            
            minZscale = std::stod(*it);
            
            if(std::next(it) != arg.cend()){it++;}else{return false;}
            if(it->find_first_not_of("+-0123456789.eE") != std::string::npos)
                return false;
            
            maxZscale = std::stod(*it);
            
            if(minZscale >= maxZscale)
            {
                std::cerr<<"[ERROR] UOption::UserParameters: Invalid zscale range ["<<minZscale<<","<<maxZscale<<"]"<<std::endl;
                return false;
            }
            
            std::cout<<"   \033[34m• Set zscale range to:\033[0m["<<minZscale<<","<<maxZscale<<"]"<<std::endl;
        }
        else if(!it->compare("--nodraw"))
        {
            UOption::ddraw = UOption::drawing_mode::OFF;
        }
        else if(!it->compare("--roi"))
        {
            if(std::next(it) != arg.cend()){it++;}else{return false;}
            ROI.insert(std::pair<char,std::pair<size_t,size_t> >('x',std::pair<size_t,size_t>(0,0)));
            ROI.find('x')->second.first = static_cast<size_t>(std::stoi(*it));
            
            if(std::next(it) != arg.cend()){it++;}else{return false;}
            ROI.find('x')->second.second = static_cast<size_t>(std::stoi(*it));
            
            if(std::next(it) != arg.cend()){it++;}else{return false;}
            ROI.insert(std::pair<char,std::pair<size_t,size_t> >('y',std::pair<size_t,size_t>(0,0)));
            ROI.find('y')->second.first = static_cast<size_t>(std::stoi(*it));
            
            if(std::next(it) != arg.cend()){it++;}else{return false;}
            ROI.find('y')->second.second = static_cast<size_t>(std::stoi(*it));

            std::cout<<"   \033[34m• Restrict to ROI:\033[0m["<<ROI['x'].first<<","<<ROI['x'].second<<"]["
                <<ROI['y'].first<<","<<ROI['y'].second<<"]"<<std::endl;
            
        }
        else if(!it->compare("-print") || !it->compare("--print") || !it->compare("+print") || !it->compare("++print"))
        {
            do_Print = true;
        }
        else
            return false;
        
        return true;
    }
    
    inline const size_t&       Rebin() const {return RBIN;}
    inline const rmode&        Mode()  const {return MODE;}
    inline const int16_t&      Header()const {return HDU;}
    inline bool doPrint() const {return do_Print;}
    inline bool ArrayInfo() const {return a_info;}
    inline bool HduExplicit() const {return hdu_explicit;}

};

UOption::drawing_mode UOption::ddraw = UOption::drawing_mode::LINEAR;

#ifdef __WITH_QT__

#ifdef __WITH_ASDF__

// enumerateHDUs(FITSmanager&) for a file of either format. Unlike a FITS file, the first HDU of an
// ASDF file can be a table, so XTENSION is what decides the type here, whatever the index.
static std::vector<HduInfo> enumerateHDUs(DataFile& file)
{
    if (FITSmanager* fits = file.fitsFile())
        return enumerateHDUs(*fits);

    std::vector<HduInfo> result;
    const int n = file.NumberOfHeader();
    for (int i = 1; i <= n; ++i) {
        try {
            const std::shared_ptr<FITShdu> hdr = file.GetHeaderAtIndex(i);

            HduInfo info;
            info.index   = i;
            info.extname = hdr->Exists("EXTNAME") ? hdr->GetValueForKey("EXTNAME") : "HDU" + std::to_string(i);
            info.type    = HduInfo::Type::IMG;
            if (hdr->Exists("XTENSION") && hdr->GetValueForKey("XTENSION").find("BINTABLE") != std::string::npos)
                info.type = HduInfo::Type::BTBL;
            info.naxis   = hdr->Exists("NAXIS")  ? (int)hdr->GetUInt64ValueForKey("NAXIS")  : 0;
            info.naxis3  = hdr->Exists("NAXIS3") ? (int)hdr->GetUInt64ValueForKey("NAXIS3") : 0;

            result.push_back(info);
        }
        catch (const FITSexception& e) {
            std::cerr << "[WARNING] enumerateHDUs: HDU " << i << ": " << e.what() << std::endl;
        }
    }
    return result;
}

#endif

std::pair<double,double> getYDataRange(const QCPGraph *plt)
{
    
    const QCPRange visibleX = plt->keyAxis()->range();
    QSharedPointer<QCPGraphDataContainer> data = plt->data();

    auto itLow  = data->findBegin(visibleX.lower, true);
    auto itHigh = data->findEnd(visibleX.upper, true);

    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();

    for (auto it = itLow; it != itHigh; ++it)
    {
        const double v = it->value;
        if (!std::isfinite(v)) continue;
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }

    const double span = vmax - vmin;
    const double pad = std::max(1e-9, span * 0.05);
    
    return std::pair<double,double>(vmin - pad,vmax + pad);
}

void mouseHasMove(QCustomPlot *plt,const std::shared_ptr<const FITScube>& in, QCPGraph *x,QCPGraph *y)
{
    // Each profile below replots the whole window: while a drag pans the image, that is two more
    // full redraws per mouse move for a pixel the cursor doesn't stay on.
    if (QGuiApplication::mouseButtons() != Qt::NoButton)
        return;

    QCPAxisRect* img = dynamic_cast<QCPAxisRect*>(plt->plotLayout()->element(0,0));
    if (!img) return;

    // Only refresh when mouse is inside the image axis rect
    const QPoint cursor = plt->mapFromGlobal(QCursor::pos());
    QCPAxisRect* hit = plt->axisRectAt(cursor);
    if (hit != img)
        return;

    const QCPRange xrange = img->axis(QCPAxis::atBottom)->range();
    const QCPRange yrange = img->axis(QCPAxis::atLeft)->range();
    std::vector<double> px_coo({img->axis(QCPAxis::atBottom)->pixelToCoord(cursor.x()),img->axis(QCPAxis::atLeft)->pixelToCoord(cursor.y())});

    if (px_coo[0] < 0 || px_coo[0] >= in->Size(1) || px_coo[0]< xrange.lower || px_coo[0]> xrange.upper ||
        px_coo[1] < 0 || px_coo[1] >= in->Size(2) || px_coo[1]< yrange.lower || px_coo[1]> yrange.upper) 
        return;

    QString valStr = QString("N/A");

    try
    {
        QtFITSviewer::QtProfileX(plt, in,static_cast<size_t>(px_coo[1]+0.5),x);
        QtFITSviewer::QtProfileY(plt, in,static_cast<size_t>(px_coo[0]+0.5),y, true);

        std::pair<double,double> range_x = getYDataRange(x);
        std::pair<double,double> range_y = getYDataRange(y);

        const QCPRange current_x = x->valueAxis()->range();
        const QCPRange target_x(range_x.first, range_x.second);

        const QCPRange current_y = y->valueAxis()->range();
        const QCPRange target_y(range_y.first, range_y.second);

        // Avoid needless replots if unchanged
        if (std::abs(current_x.lower - target_x.lower) > 1e-9 || std::abs(current_x.upper - target_x.upper) > 1e-9)
                x->valueAxis()->setRange(target_x);

        // Avoid needless replots if unchanged
        if (std::abs(current_y.lower - target_y.lower) > 1e-9 || std::abs(current_y.upper - target_y.upper) > 1e-9)
                y->valueAxis()->setRange(target_y);

        // --- NEW: update main window status bar with pixel coords and pixel value ---
        // compute integer pixel indices (col=X, row=Y)
        const size_t ix = static_cast<size_t>(px_coo[0] + 0.5);
        const size_t iy = static_cast<size_t>(px_coo[1] + 0.5);

        // Placeholder accessor: replace with your FITScube method (example: GetValueAt, At, operator(), ...)
        double pixVal = 0.0;
        bool got = false;

        pixVal = in->DoubleValueAtPixel({ix, iy});
        got = true;

        if(got)
            valStr = QString::number(pixVal, 'g', 6);

        QString wcs_msg = QString("");
        if(in->HasWorldCoordinates() && in->GetDimension() == 2)
        {
            // ix, iy are 0-based (colour-map cells), the convention FITScube works in.
            // The FITS WCS of the header, or another WCS of the image (the GWCS of an ASDF image).
            worldCoords wc = in->WorldCoordinates(std::vector<size_t>({ix, iy}));
            std::vector< std::string > ctype = in->WorldAxisNames();
            std::vector< std::string > cunit = in->WorldAxisUnits();
            wc.resize(2, std::numeric_limits<double>::quiet_NaN());
            ctype.resize(2);    // A WCS may leave its axes unnamed
            cunit.resize(2);

            std::string xstr = (ctype[0].size() > 0)?ctype[0]:"\u03b1"; 
            if(xstr.find('-') != std::string::npos) xstr.resize(xstr.find_first_of('-'));
            std::string ystr = (ctype[1].size() > 0)?ctype[1]:"\u03b4"; 
            if(ystr.find('-') != std::string::npos) ystr.resize(ystr.find_first_of('-'));

            std::string xUnit = (cunit[0].size() > 0)?" "+cunit[0]:""; if(xUnit == " deg" || xUnit == " DEG" ) xUnit = "\u00b0";
            std::string yUnit = (cunit[1].size() > 0)?" "+cunit[1]:""; if(yUnit == " deg" || yUnit == " DEG" ) yUnit = "\u00b0";

            // Degrees with 7 decimals: 0.36 mas, ~1/1000 of a NISP pixel (0.3"). 3 decimals
            // (3.6") was coarser than 10 pixels. Other units: 10 significant digits.
            auto fmtWorld = [](double v, const std::string& unit)->QString
            {
                return (unit == "°") ? QString::number(v, 'f', 7) : QString::number(v, 'g', 10);
            };

            wcs_msg = QString("%1=%2%3  %4=%5%6 ")
                          .arg(QString::fromStdString(xstr)).arg(fmtWorld(wc[0], xUnit)).arg(QString::fromStdString(xUnit))
                          .arg(QString::fromStdString(ystr)).arg(fmtWorld(wc[1], yUnit)).arg(QString::fromStdString(yUnit));
        }

        // find top-level window (the QtFITSviewer is a QMainWindow) and write to its status bar
        QWidget *top = plt->window();
        if(top)
        {
            QMainWindow *mainWin = qobject_cast<QMainWindow*>(top);
            if(mainWin && mainWin->statusBar())
            {
                QString msg = QString("[%1 ; %2] ").arg(ix).arg(iy)+wcs_msg+QString("val=%3").arg(valStr);
                mainWin->statusBar()->showMessage(msg);
            }
        }
    }
    catch(std::out_of_range& e)
    {
        if((DSL::verbose&DSL::verboseLevel::VERBOSE_DEBUG) == DSL::verboseLevel::VERBOSE_DEBUG)
        std::cerr << "[WARNING] mouseHasMove: " << e.what() << std::endl;
        valStr = QString("N/A");
    }


}

void PrintInfo(QCustomPlot *plt, const std::shared_ptr<const FITScube>& in)
{
    if(!in)
        return;

    std::string unit = "";
    if(in->HDU().Exists("BUNIT"))
        unit += in->HDU().GetValueForKey("BUNIT");
    else
        unit += " ARB. UNIT";

    double mean = in->GetMean();
    double median = in->GetMedian();
    double stddev = in->GetStdDev();
    double min   = in->GetMinimum();
    double max   = in->GetMaximum();
    double p05   = in->Getpercentil(0.05);
    double p25   = in->Getpercentil(0.25);
    double p75   = in->Getpercentil(0.75);
    double p95   = in->Getpercentil(0.95);
    
    QString info = QString::fromStdString(
         "Min     : "+DST::to_string(min,3)+" "+unit+"   "
        +"Max     : "+DST::to_string(max,3)+" "+unit+"   "
        +"Mean    : "+DST::to_string(mean,3)+" "+unit+"\n"
        +"5th %   : "+DST::to_string(p05,3)+" "+unit+"   "
        +"Median  : "+DST::to_string(median,3)+" "+unit+"   "
        +"95th %  : "+DST::to_string(p95,3)+" "+unit+"\n"
        +"Stddev  : "+DST::to_string(stddev,3)+" "+unit+"   "
        +"Kurtosis: "+DST::to_string(in->GetKurtosis(),3)+" "+unit+"   "
        +"Skewness: "+DST::to_string(in->GetSkewness(),3)+" "+unit+"\n"
        );

    QMetaObject::invokeMethod(plt, [info,p05,p25,median,p75,p95,mean,plt ]() mutable
                                            {
                                                {
                                                    std::lock_guard<std::mutex> lk(QtFITS::io_mutex);
                                                    // add the text element to the layout (keeps it in the (1,3) cell)

                                                    QCPTextElement *imgStatElement = new QCPTextElement(plt);
                                                    imgStatElement->setTextFlags(Qt::AlignLeft|Qt::TextWordWrap);
                                                    imgStatElement->setText(info);

                                                    imgStatElement->setFont(QFont(QToolTip::font().family(), 12));

                                                    QCPItemStraightLine *p05fLine = new QCPItemStraightLine(plt);
                                                    p05fLine->setPen(QPen(Qt::blue, 1, Qt::DotLine));
                                                    p05fLine->point1->setCoords(p05, 0); // point1 is on the line
                                                    p05fLine->point2->setCoords(p05, 1); //

                                                    QCPItemStraightLine *p95fLine = new QCPItemStraightLine(plt);
                                                    p95fLine->setPen(QPen(Qt::blue, 1, Qt::DotLine));
                                                    p95fLine->point1->setCoords(p95, 0); // point1 is on the line
                                                    p95fLine->point2->setCoords(p95, 1); //

                                                    QCPItemStraightLine *p25fLine = new QCPItemStraightLine(plt);
                                                    p25fLine->setPen(QPen(Qt::blue, 1, Qt::DashLine));
                                                    p25fLine->point1->setCoords(p25, 0); // point1 is on the line
                                                    p25fLine->point2->setCoords(p25, 1); //

                                                    QCPItemStraightLine *p75fLine = new QCPItemStraightLine(plt);
                                                    p75fLine->setPen(QPen(Qt::blue, 1, Qt::DashLine));
                                                    p75fLine->point1->setCoords(p75, 0); // point1 is on the line
                                                    p75fLine->point2->setCoords(p75, 1); //

                                                    QCPItemStraightLine *p50fLine = new QCPItemStraightLine(plt);
                                                    p50fLine->setPen(QPen(Qt::blue, 2));
                                                    p50fLine->point1->setCoords(median, 0); // point1 is on the line
                                                    p50fLine->point2->setCoords(median, 1); //

                                                    QCPItemStraightLine *meanfLine = new QCPItemStraightLine(plt);
                                                    meanfLine->setPen(QPen(Qt::red, 2));
                                                    meanfLine->point1->setCoords(mean, 0); // point1 is on the line
                                                    meanfLine->point2->setCoords(mean, 1); //

                                                    plt->plotLayout()->addElement(1, 3, imgStatElement);
                                                    QCPLayoutElement *ell = plt->plotLayout()->element(1, 3);
                                                    ell->setMinimumSize(100, 50);
                                                    //ell->setMaximumSize(1200, 1200);
                                                    

                                                    // Try to locate the histogram axis rect we created at layout (0,3)
                                                    QCPAxisRect *histoRect = nullptr;
                                                    QCPLayoutElement *el = plt->plotLayout()->element(0, 3);
                                            
                                                    if (el) {
                                                        // element at (0,3) is a QCPLayoutGrid (histo_layout). The histo axis rect
                                                        // was added as grid->element(0,0) in Draw().
                                                        QCPLayoutGrid *grid = dynamic_cast<QCPLayoutGrid*>(el);
                                                        if (grid) {
                                                            QCPLayoutElement *child = grid->element(0, 0);
                                                            histoRect = dynamic_cast<QCPAxisRect*>(child);
                                                        } else {
                                                            // sometimes the element may be directly the axis rect
                                                            histoRect = dynamic_cast<QCPAxisRect*>(el);
                                                        }
                                                    }

                                                    if (histoRect) {
                                                        // bind the item's positions to the histo axes so the line appears inside that panel
                                                        auto xAxis = histoRect->axis(QCPAxis::atBottom);
                                                        auto yAxis = histoRect->axis(QCPAxis::atLeft);
                                                        if (xAxis && yAxis) {
                                                            p05fLine->point1->setAxes(xAxis, yAxis);
                                                            p25fLine->point1->setAxes(xAxis, yAxis);
                                                            p50fLine->point1->setAxes(xAxis, yAxis);
                                                            p75fLine->point1->setAxes(xAxis, yAxis);
                                                            p95fLine->point1->setAxes(xAxis, yAxis);
                                                            meanfLine->point1->setAxes(xAxis, yAxis);
                                                            p05fLine->point2->setAxes(xAxis, yAxis);
                                                            p25fLine->point2->setAxes(xAxis, yAxis);
                                                            p50fLine->point2->setAxes(xAxis, yAxis);
                                                            p75fLine->point2->setAxes(xAxis, yAxis);
                                                            p95fLine->point2->setAxes(xAxis, yAxis);
                                                            meanfLine->point2->setAxes(xAxis, yAxis);
                                                        }
                                                        // clip the item to that axis rect so it doesn't draw outside the panel
                                                        p05fLine->setClipToAxisRect(true);
                                                        p25fLine->setClipToAxisRect(true);
                                                        p50fLine->setClipToAxisRect(true);
                                                        p75fLine->setClipToAxisRect(true);
                                                        p95fLine->setClipToAxisRect(true);
                                                        meanfLine->setClipToAxisRect(true);
                                                        p05fLine->setClipAxisRect(histoRect);
                                                        p25fLine->setClipAxisRect(histoRect);
                                                        p50fLine->setClipAxisRect(histoRect);
                                                        p75fLine->setClipAxisRect(histoRect);
                                                        p95fLine->setClipAxisRect(histoRect);
                                                        meanfLine->setClipAxisRect(histoRect);
                                                    } else {
                                                        // fallback: leave previous behaviour (global coords)
                                                        p05fLine->setClipToAxisRect(false);
                                                        p25fLine->setClipToAxisRect(false);
                                                        p50fLine->setClipToAxisRect(false);
                                                        p75fLine->setClipToAxisRect(false);
                                                        p95fLine->setClipToAxisRect(false);
                                                        meanfLine->setClipToAxisRect(false);
                                                        
                                                    }

                                                    // item is parented to plt (created with plt) so no need to add to layout.
                                                    plt->replot();
                                                } // io_mutex released here
                                            
                                                {
                                                    std::lock_guard<std::mutex> dk(QtFITS::msg_m);
                                                    if (!QtFITS::finished)
                                                         QtFITS::finished |= true;
                                                    QtFITS::msg_cv.notify_one(); // notify while holding msg_m (or move notify outside if preferred)
                                                    QtFITS::finished &= false;
                                                } // msg_m released here
                                            }
                                        );
    
}

void onResize(QCustomPlot *plt)
{
    QMetaObject::invokeMethod(plt, [plt ]() mutable
    {
        {
            std::lock_guard<std::mutex> lk(QtFITS::io_mutex);

            auto img =  dynamic_cast<QCPAxisRect*>(plt->plotLayout()->element(0,0));
            QCPAxis *x = img->axis(QCPAxis::atBottom);
            QCPAxis *y = img->axis(QCPAxis::atLeft);
            
            //const int width = plt->window()->width();
            //const int height = plt->window()->height();

            //img->setMinimumSize(height * 2/3, height * 2/3);

            //auto xproj = plt->plotLayout()->element(1,0);
            //xproj->setMinimumSize(img->width(), height *1/3);

            //auto yproj = plt->plotLayout()->element(0,2);
            //yproj->setMinimumSize(height *1/3, img->height());

            //auto hist = plt->plotLayout()->element(0,3);
            //hist->setMinimumSize(width - (height), img->height());

            y->setScaleRatio(x, 1.0);
            plt->replot();
        }

        {
            std::lock_guard<std::mutex> dk(QtFITS::msg_m);
            if (!QtFITS::finished)
                 QtFITS::finished |= true;
            QtFITS::msg_cv.notify_one(); // notify while holding msg_m (or move notify outside if preferred)
            QtFITS::finished &= false;
        } // msg_m released here
    });
}

void DrawWCSgrid(QCustomPlot *plt, const std::shared_ptr<const FITScube>& in)
{
    if(!in)
        return;

    if(!in->HasWorldCoordinates())
        return;

    // A FITS WCS that is only the pixel grid has no sky to draw; checked only when it is the WCS in use (an
    // ASDF image with a GWCS names its axes after the GWCS, whatever its header holds).
    if( in->getNumberOfWCS() > 0 && in->WorldAxisNames() == in->getWCS().CTYPE(0) &&
        in->getWCS().CRPIX(0,1) == 0 && in->getWCS().CRPIX(0,2) == 0 &&
        in->getWCS().CRVAL(0,1) == 0 && in->getWCS().CRVAL(0,2) == 0 &&
        in->getWCS().CDELT(0,1) == 1 && in->getWCS().CDELT(0,2) == 1)
        return;

    const size_t nx = in->Size(1);
    const size_t ny = in->Size(2);

    // Direct grid-line tracing — no full-resolution WCS map, no contours.
    //
    // Approach: for each target coordinate value, scan columns (for Dec lines)
    // or rows (for RA lines) and interpolate the pixel position where the WCS
    // coordinate crosses the target.  Total WCS calls ≈ n_lines × image_dim/stride
    // × N_PROBE ≈ a few thousand, vs O(nx·ny) for the old contour approach.
    //
    // The WCS cache avoids redundant evaluations when multiple grid lines probe
    // the same (column, probe_y) or (probe_x, row) locations.

    static constexpr int    N_PROBE = 8;   // probe points per col/row to bracket crossings
    static constexpr size_t N_MAX   = 256; // max columns/rows to trace per axis direction

    const size_t col_stride = std::max(size_t(1), (nx + N_MAX - 1) / N_MAX);
    const size_t row_stride = std::max(size_t(1), (ny + N_MAX - 1) / N_MAX);

    // WCS cache keyed on the 1-D pixel index so repeated probes come for free.
    // std::map is used (no extra include needed) since n is small (≤ N_MAX×N_PROBE×2).
    std::map<size_t, worldCoords> cache;

    auto wcs = [&](size_t ix, size_t iy) -> const worldCoords& {
        size_t key = ix + iy * nx;
        auto it = cache.find(key);
        if (it == cache.end())
            it = cache.emplace(key, in->WorldCoordinates(key)).first;
        return it->second;
    };

    // Sample corners + centre to determine the visible RA/Dec extent.
    const auto& w00  = wcs(0,    0);
    const auto& w10  = wcs(nx-1, 0);
    const auto& w01  = wcs(0,    ny-1);
    const auto& w11  = wcs(nx-1, ny-1);
    const auto& wMid = wcs(nx/2, ny/2);

    double ra_min  = std::min({w00[0], w10[0], w01[0], w11[0], wMid[0]});
    double ra_max  = std::max({w00[0], w10[0], w01[0], w11[0], wMid[0]});
    double dec_min = std::min({w00[1], w10[1], w01[1], w11[1], wMid[1]});
    double dec_max = std::max({w00[1], w10[1], w01[1], w11[1], wMid[1]});

    // Pick a "nice" grid step that gives 4–6 lines per axis.
    auto niceStep = [](double span) -> double {
        if (span <= 0.0) return 1.0;
        double raw = span / 5.0;
        double e   = std::pow(10.0, std::floor(std::log10(raw)));
        double m   = raw / e;
        if (m <= 1.5) return e;
        if (m <= 3.5) return 2.0 * e;
        if (m <= 7.5) return 5.0 * e;
        return 10.0 * e;
    };

    const double ra_step  = niceStep(ra_max  - ra_min);
    const double dec_step = niceStep(dec_max - dec_min);

    // Format a coordinate value with enough decimal places for the step size.
    auto fmtCoord = [](double val, double step) -> std::string {
        int dp = std::max(0, (int)std::ceil(-std::log10(std::max(step, 1e-10))));
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(dp) << val << "\xc2\xb0"; // UTF-8 °
        return oss.str();
    };

    if(_cancel_wcs) return;

    // Each grid line stores its polyline plus the position of its label (midpoint
    // of the first valid segment) and the label text.
    struct GridLine {
        std::vector<double> xs, ys;
        double lx = std::numeric_limits<double>::quiet_NaN();
        double ly = std::numeric_limits<double>::quiet_NaN();
        std::string label;
        size_t npts = 0; // valid point count, used to pick the midpoint
    };
    std::vector<GridLine> grid_lines;

    const double nan_ = std::numeric_limits<double>::quiet_NaN();

    // ---- Lines of constant Dec (scan columns, find y-crossing) ----
    {
        double d = std::ceil((dec_min - 1e-9 * dec_step) / dec_step) * dec_step;
        for (; d <= dec_max + dec_step * 1e-6; d += dec_step)
        {
            if(_cancel_wcs) return;
            GridLine gl;
            gl.label = fmtCoord(d, dec_step);
            bool hadPoint = false;

            for (size_t ix = 0; ix < nx; ix += col_stride)
            {
                double found_y = nan_;
                double prev_dec = nan_;
                double prev_fy  = 0.0;

                for (int k = 0; k <= N_PROBE; ++k)
                {
                    double fy  = double(k) / double(N_PROBE) * double(ny - 1);
                    size_t iy  = static_cast<size_t>(std::round(fy));
                    double dec = wcs(ix, iy)[1];

                    if (k > 0 && std::isfinite(prev_dec) && std::isfinite(dec) &&
                        (dec - prev_dec) != 0.0 &&
                        (prev_dec - d) * (dec - d) <= 0.0)
                    {
                        double t = (d - prev_dec) / (dec - prev_dec);
                        found_y = prev_fy + t * (fy - prev_fy);
                        break;
                    }
                    prev_dec = dec;
                    prev_fy  = fy;
                }

                if (std::isfinite(found_y))
                {
                    double px = static_cast<double>(ix);
                    gl.xs.push_back(px);
                    gl.ys.push_back(found_y);
                    ++gl.npts;
                    // Track midpoint of the first contiguous segment for the label.
                    // We use every other valid point as the candidate so that the
                    // label ends up near the middle of the line.
                    if (gl.npts % 2 == 1) { gl.lx = px; gl.ly = found_y; }
                    hadPoint = true;
                }
                else if (hadPoint)
                {
                    gl.xs.push_back(nan_);
                    gl.ys.push_back(nan_);
                    hadPoint = false;
                }
            }
            if (gl.xs.size() > 1)
                grid_lines.push_back(std::move(gl));
        }
    }

    // ---- Lines of constant RA (scan rows, find x-crossing) ----
    {
        double r = std::ceil((ra_min - 1e-9 * ra_step) / ra_step) * ra_step;
        for (; r <= ra_max + ra_step * 1e-6; r += ra_step)
        {
            if(_cancel_wcs) return;
            GridLine gl;
            gl.label = fmtCoord(r, ra_step);
            bool hadPoint = false;

            for (size_t iy = 0; iy < ny; iy += row_stride)
            {
                double found_x = nan_;
                double prev_ra = nan_;
                double prev_fx = 0.0;

                for (int k = 0; k <= N_PROBE; ++k)
                {
                    double fx = double(k) / double(N_PROBE) * double(nx - 1);
                    size_t ix = static_cast<size_t>(std::round(fx));
                    double ra = wcs(ix, iy)[0];

                    if (k > 0 && std::isfinite(prev_ra) && std::isfinite(ra) &&
                        (ra - prev_ra) != 0.0 &&
                        (prev_ra - r) * (ra - r) <= 0.0)
                    {
                        double t = (r - prev_ra) / (ra - prev_ra);
                        found_x = prev_fx + t * (fx - prev_fx);
                        break;
                    }
                    prev_ra = ra;
                    prev_fx = fx;
                }

                if (std::isfinite(found_x))
                {
                    double py = static_cast<double>(iy);
                    gl.xs.push_back(found_x);
                    gl.ys.push_back(py);
                    ++gl.npts;
                    if (gl.npts % 2 == 1) { gl.lx = found_x; gl.ly = py; }
                    hadPoint = true;
                }
                else if (hadPoint)
                {
                    gl.xs.push_back(nan_);
                    gl.ys.push_back(nan_);
                    hadPoint = false;
                }
            }
            if (gl.xs.size() > 1)
                grid_lines.push_back(std::move(gl));
        }
    }

    if(_cancel_wcs || grid_lines.empty()) return;

    QMetaObject::invokeMethod(plt, [plt, gls = std::move(grid_lines)]() mutable
    {
        {
            std::lock_guard<std::mutex> lk(QtFITS::io_mutex);
            for (const auto& gl : gls)
            {
                // Draw the polyline
                if (!gl.xs.empty())
                {
                    auto* curve = new QCPCurve(plt->xAxis, plt->yAxis);
                    curve->setPen(QPen(Qt::green, 1.0));
                    curve->setScatterStyle(QCPScatterStyle::ssNone);
                    curve->setLayer("wcs_grid");
                    for (size_t i = 0; i < gl.xs.size(); ++i)
                        curve->addData(gl.xs[i], gl.ys[i]);
                }

                // Add a coordinate label at the midpoint of the line
                if (std::isfinite(gl.lx) && std::isfinite(gl.ly) && !gl.label.empty())
                {
                    auto* lbl = new QCPItemText(plt);
                    lbl->setPositionAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
                    lbl->position->setType(QCPItemPosition::ptPlotCoords);
                    lbl->position->setCoords(gl.lx, gl.ly);
                    lbl->setText(QString::fromStdString(gl.label));
                    lbl->setColor(Qt::green);
                    QFont lblFont = QFontDatabase::systemFont(QFontDatabase::FixedFont);
                    lblFont.setPointSize(7);
                    lbl->setFont(lblFont);
                    lbl->setPen(Qt::NoPen);
                    // Semi-transparent dark background so the text is readable
                    // over the image regardless of pixel intensity.
                    lbl->setBrush(QBrush(QColor(0, 0, 0, 140)));
                    lbl->setPadding(QMargins(2, 1, 2, 1));
                    lbl->setLayer("wcs_grid");
                }
            }
            plt->replot();
        }
        {
            std::lock_guard<std::mutex> dk(QtFITS::msg_m);
            if (!QtFITS::finished) QtFITS::finished |= true;
            QtFITS::msg_cv.notify_one();
            QtFITS::finished &= false;
        }
    });
}

QCustomPlot* Draw(const std::shared_ptr<const FITScube>& in)
{
    if(UOption::ddraw == UOption::drawing_mode::OFF)
        return NULL;
    
    QtFITSviewer *mw = new QtFITSviewer(400, 250, 1200, 900);
    mw->setWindowTitle(in->GetName().c_str());
    mw->linkResizeEvent(onResize);

    //size_t xsize = std::min(static_cast<size_t>(900),in->Size(1)+200);
    //size_t ysize = std::min(static_cast<size_t>(900),in->Size(2)+200);

    QCustomPlot *customPlot = mw->GetPlot();
    
    QtFITSviewer::QtDisplay(customPlot,in);
    
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes );
    QCPAxisRect *image = customPlot->axisRect(0);
    image->setBackgroundScaledMode(Qt::AspectRatioMode::KeepAspectRatioByExpanding);
    
    QCPLayoutGrid *xproj_layout = new QCPLayoutGrid;
    customPlot->plotLayout()->addElement(1, 0, xproj_layout); // row 1, col 0: x projection

    QCPLayoutGrid *yproj_layout = new QCPLayoutGrid;
    customPlot->plotLayout()->addElement(0, 2, yproj_layout); // row 0, col 2: y projection (col 1 is color scale)

    QCPLayoutGrid *histo_layout = new QCPLayoutGrid;
    customPlot->plotLayout()->addElement(0, 3, histo_layout); // row 0, col 3: histogram

    // Col 0=image(6), col 1=colorscale(narrow), col 2=yproj(1), col 3=histo(1)
    customPlot->plotLayout()->setColumnStretchFactor(0, 6);
    customPlot->plotLayout()->setColumnStretchFactor(1, 0.3);
    customPlot->plotLayout()->setColumnStretchFactor(2, 1);
    customPlot->plotLayout()->setColumnStretchFactor(3, 3);

    // Row 0=image(3), row 1=xproj(1)
    customPlot->plotLayout()->setRowStretchFactor(0, 3);
    customPlot->plotLayout()->setRowStretchFactor(1, 1);

    // prepare axis rects that will be placed in the sublayout:
    QCPAxisRect *yproj = new QCPAxisRect(customPlot, false);
    QCPAxisRect *xproj = new QCPAxisRect(customPlot, false);
    QCPAxisRect *histo = new QCPAxisRect(customPlot, false);

    image->setMinimumSize(200, 200);
    //image->setMaximumSize(1900, 1900);
    yproj->setMinimumSize(150, 200);
    xproj->setMinimumSize(200, 150);
    //xproj->setMaximumSize(QWIDGETSIZE_MAX, 100);
    histo->setMinimumSize(200, 200);
    
    yproj_layout->addElement(0, 0, yproj);
    xproj_layout->addElement(0, 0, xproj);
    histo_layout->addElement(0, 0, histo);
    
    
    xproj->addAxes(QCPAxis::atBottom | QCPAxis::atLeft);
    xproj->axis(QCPAxis::atBottom)->setRange(image->axis(QCPAxis::atBottom)->range());
    xproj->axis(QCPAxis::atLeft)->ticker()->setTickCount(2);
    xproj->axis(QCPAxis::atBottom)->ticker()->setTickCount(2);
    xproj->axis(QCPAxis::atLeft)->setRange(in->GetMinimum(),in->GetMaximum());
    
    yproj->addAxes(QCPAxis::atBottom | QCPAxis::atRight);
    yproj->axis(QCPAxis::atRight)->setRange(image->axis(QCPAxis::atLeft)->range());
    yproj->axis(QCPAxis::atRight)->ticker()->setTickCount(2);
    yproj->axis(QCPAxis::atBottom)->ticker()->setTickCount(2);
    yproj->axis(QCPAxis::atBottom)->setRangeReversed(true);
    yproj->axis(QCPAxis::atBottom)->setRange(in->GetMinimum(),in->GetMaximum());
    
    histo->addAxes(QCPAxis::atBottom | QCPAxis::atRight | QCPAxis::atLeft);
    histo->axis(QCPAxis::atBottom)->setRange(in->GetMinimum()*(1-0.05), in->GetMaximum()*(1+0.05));
    histo->axis(QCPAxis::atBottom)->setScaleType(QCPAxis::stLinear);
    histo->axis(QCPAxis::atBottom)->ticker()->setTickCount(2);
    histo->axis(QCPAxis::atRight)->setRange(-0.05, 1.05);
    histo->axis(QCPAxis::atRight)->setScaleType(QCPAxis::stLinear);
    histo->axis(QCPAxis::atRight)->ticker()->setTickCount(2);
    histo->axis(QCPAxis::atLeft)->setScaleType(QCPAxis::stLogarithmic);
    
    histo->setRangeDragAxes({histo->axis(QCPAxis::atBottom),histo->axis(QCPAxis::atRight),histo->axis(QCPAxis::atLeft)});
    histo->setRangeZoomAxes({histo->axis(QCPAxis::atBottom),histo->axis(QCPAxis::atRight),histo->axis(QCPAxis::atLeft)});
    
    
    QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
    histo->axis(QCPAxis::atLeft)->setTicker(logTicker);
    
    foreach (QCPAxis *axis, yproj->axes())
        axis->setLayer("axes");
    
    foreach (QCPAxis *axis, xproj->axes())
        axis->setLayer("axes");
    
    foreach (QCPAxis *axis, histo->axes())
        axis->setLayer("axes");
    
    QCPGraph *x_projection = customPlot->addGraph(xproj->axis(QCPAxis::atBottom), xproj->axis(QCPAxis::atLeft));
    x_projection->setLineStyle(QCPGraph::lsStepCenter);
    x_projection->setPen(QPen(QColor("#FFA100"), 1.5));
    
    QCPGraph *y_projection = customPlot->addGraph(yproj->axis(QCPAxis::atRight), yproj->axis(QCPAxis::atBottom));
    y_projection->setLineStyle(QCPGraph::lsStepCenter);
    y_projection->setPen(QPen(QColor("#FFA100"), 1.5));
    
    QCPGraph *px_pdf = customPlot->addGraph(histo->axis(QCPAxis::atBottom), histo->axis(QCPAxis::atLeft));
    px_pdf->setLineStyle(QCPGraph::lsStepCenter);
    px_pdf->setPen(QPen(QColor(255, 100, 0)));

    QCPGraph *px_cdf = customPlot->addGraph(histo->axis(QCPAxis::atBottom), histo->axis(QCPAxis::atRight));
    px_cdf->setLineStyle(QCPGraph::lsLine);
    px_cdf->setPen(QPen(QColor(30, 40, 255, 150)));
    //customPlot->replot();


    if(minZscale == minZscale && maxZscale  == maxZscale)
    {
        QCPColorScale *colorScale = NULL;
        if( customPlot->plotLayout()->hasElement (0, 1) )
        {
            colorScale = dynamic_cast<QCPColorScale*>(customPlot->plotLayout()->element(0,1));
            colorScale->setDataRange((QCPRange(minZscale,maxZscale)));
        }
    }
    

    _draw_histo= std::thread([=]() {QtFITSviewer::QtHistogram(customPlot, in, static_cast<size_t>(histo_bins), px_pdf, px_cdf);});
    _print_info = std::thread([=]() {PrintInfo(customPlot,in);});
    

    customPlot->connect(image->axis(QCPAxis::atBottom), SIGNAL(rangeChanged(QCPRange)), xproj->axis(QCPAxis::atBottom), SLOT(setRange(QCPRange)));
    
    customPlot->connect(image->axis(QCPAxis::atLeft), SIGNAL(rangeChanged(QCPRange)), yproj->axis(QCPAxis::atRight), SLOT(setRange(QCPRange)));
    
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes | QCP::iSelectOther);
    customPlot->setNoAntialiasingOnDrag(true);
    customPlot->setPlottingHint(QCP::phImmediateRefresh, false);
    //customPlot->replot();
    
    auto f = std::bind(mouseHasMove,customPlot,in,x_projection,y_projection);
    QAbstractAnimation::connect(customPlot, &QCustomPlot::mouseMove, f);

    image->axis(QCPAxis::atLeft)->setScaleRatio(image->axis(QCPAxis::atBottom), 1.0);
    QObject::connect(customPlot, &QCustomPlot::beforeReplot, [image]() {
        image->axis(QCPAxis::atLeft)->setScaleRatio(image->axis(QCPAxis::atBottom), 1.0);
    });

    customPlot->addLayer("wcs_grid", customPlot->layer("main"), QCustomPlot::limAbove);

    mw->SetPlot(customPlot);
    mw->show();

    _draw_grid= std::thread([=]() {DrawWCSgrid(customPlot, in);});

    return customPlot;
}

// ---------------------------------------------------------------------------
// ArrowSliderStyle — downward-pointing triangle handle with tick marks
// ---------------------------------------------------------------------------

class ArrowSliderStyle : public QProxyStyle
{
public:
    explicit ArrowSliderStyle(QStyle* base = nullptr) : QProxyStyle(base) {}

    void drawComplexControl(ComplexControl cc, const QStyleOptionComplex* opt,
                             QPainter* p, const QWidget* w) const override
    {
        if (cc != CC_Slider) {
            QProxyStyle::drawComplexControl(cc, opt, p, w);
            return;
        }
        const auto* so = qstyleoption_cast<const QStyleOptionSlider*>(opt);
        if (!so) { QProxyStyle::drawComplexControl(cc, opt, p, w); return; }

        // Draw groove and tick marks via base style (skip the default handle)
        {
            QStyleOptionSlider tmp = *so;
            tmp.subControls &= ~QStyle::SC_SliderHandle;
            QProxyStyle::drawComplexControl(cc, &tmp, p, w);
        }

        // Draw downward-pointing triangle handle
        if (so->subControls & SC_SliderHandle) {
            const QRect hr = subControlRect(cc, opt, SC_SliderHandle, w);
            const int cx  = hr.center().x();
            const int top = hr.top();
            const int bot = hr.bottom();
            const int hw  = std::min(hr.width() / 2, 9);

            QPolygon tri;
            tri << QPoint(cx - hw, top)
                << QPoint(cx + hw, top)
                << QPoint(cx,      bot);

            const bool active = (so->state & QStyle::State_MouseOver) ||
                                (so->state & QStyle::State_Sunken);

            p->save();
            p->setRenderHint(QPainter::Antialiasing, true);
            p->setPen(QPen(active ? QColor(30, 100, 220) : QColor(50, 50, 50), 1));
            p->setBrush(active ? QColor(70, 140, 255) : QColor(90, 90, 90));
            p->drawPolygon(tri);
            p->restore();
        }
    }
};

// ---------------------------------------------------------------------------
// Viewer helpers
// ---------------------------------------------------------------------------

static QWidget* notYetImpl(const std::string& desc, const HduInfo& hdu)
{
    auto* w = new QWidget;
    w->setWindowTitle(QString("HDU %1 - placeholder").arg(hdu.index));
    auto* lbl = new QLabel(
        QString("Not yet implemented.\n%1\n(HDU index: %2)")
            .arg(QString::fromStdString(desc))
            .arg(hdu.index),
        w);
    lbl->setAlignment(Qt::AlignCenter);
    lbl->setWordWrap(true);
    auto* lay = new QVBoxLayout(w);
    lay->addWidget(lbl);
    w->resize(400, 200);
    return w;
}

// ---------------------------------------------------------------------------
// Command bar helpers
// ---------------------------------------------------------------------------

// Closes `popup` when the watched widget receives QEvent::Close.
// Parented to `popup` so it is automatically uninstalled when popup is destroyed.
class PopupCloseGuard : public QObject
{
public:
    PopupCloseGuard(QWidget* popup, QObject* parent)
        : QObject(parent), m_popup(popup) {}
protected:
    bool eventFilter(QObject*, QEvent* e) override
    {
        if (e->type() == QEvent::Close)
            m_popup->close();
        return false; // do not consume the event
    }
private:
    QWidget* m_popup;
};

// Builds and returns a populated, read-only FITS header QTableWidget.
// Caller is responsible for ownership (pass a parent or add to a layout).
static QTableWidget* makeFitsHeaderTable(const FITShdu& hdu, QWidget* parent)
{
    auto* table = new QTableWidget(parent);
    table->setColumnCount(4);
    table->setHorizontalHeaderLabels({"Keyword", "Value", "Type", "Comment"});
    table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table->setWordWrap(true);
    table->verticalHeader()->setVisible(false);
    table->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
    table->horizontalHeader()->setStretchLastSection(true);

    int row = 0;
    for (auto it = hdu.begin(); it != hdu.end(); ++it, ++row)
    {
        table->insertRow(row);
        table->setItem(row, 0, new QTableWidgetItem(QString::fromStdString(it->first)));
        table->setItem(row, 1, new QTableWidgetItem(QString::fromStdString(it->second.value())));
        table->setItem(row, 2, new QTableWidgetItem(QString::fromStdString(FITSkeyword::GetDataType(it->second.type()))));
        table->setItem(row, 3, new QTableWidgetItem(QString::fromStdString(it->second.comment())));
    }
    table->resizeColumnsToContents();
    table->setColumnWidth(1, qMin(table->columnWidth(1), table->fontMetrics().averageCharWidth() * 50));
    table->resizeRowsToContents();
    return table;
}

// Opens a non-modal window listing the FITS keywords of a given HDU.
// The widget is Qt-parented to `parent` so it is destroyed when the viewer closes.
static void showFitsHeaderWidget(const FITShdu& hdu, QWidget* parent)
{
    auto* w = new QWidget(parent, Qt::Window);
    w->setWindowTitle(QString("FITS Header — %1").arg(parent ? parent->windowTitle() : ""));
    w->setAttribute(Qt::WA_DeleteOnClose);

    auto* lay = new QVBoxLayout(w);
    lay->addWidget(makeFitsHeaderTable(hdu, w));

    w->resize(700, 400);
    if (parent)
        parent->installEventFilter(new PopupCloseGuard(w, w));
    w->show();
}

// Attaches a thin command bar (QToolBar) to any QMainWindow.
// To add more buttons in the future: create a QPushButton, connect it, call toolbar->addWidget().
static void addCommandBar(QMainWindow* win, const FITShdu& hdu, QCustomPlot* customPlot = nullptr)
{
    auto* toolbar = new QToolBar("Options", win);
    toolbar->setMovable(false);
    toolbar->setFloatable(false);
    toolbar->setIconSize(QSize(16, 16));

    // --- Header button ---
    auto* headerBtn = new QPushButton("Header", toolbar);
    headerBtn->setFlat(true);
    QObject::connect(headerBtn, &QPushButton::clicked, win, [hdu, win]() {
        showFitsHeaderWidget(hdu, win);
    });
    toolbar->addWidget(headerBtn);

    // --- WCS Grid checkbox (only present when the wcs_grid layer was created) ---
    if (customPlot && customPlot->layer("wcs_grid")) {
        toolbar->addSeparator();
        auto* wcsCheck = new QCheckBox("WCS Grid", toolbar);
        wcsCheck->setChecked(true);
        QObject::connect(wcsCheck, &QCheckBox::toggled, win, [customPlot](bool checked) {
            auto* layer = customPlot->layer("wcs_grid");
            if (layer) {
                layer->setVisible(checked);
                customPlot->replot();
            }
        });
        toolbar->addWidget(wcsCheck);
    }

    // --- Save button (all image viewers) ---
    if (customPlot) {
        toolbar->addSeparator();
        auto* saveBtn = new QPushButton("Save", toolbar);
        saveBtn->setFlat(true);
        QObject::connect(saveBtn, &QPushButton::clicked, win, [customPlot, win]() {
            QString path = QFileDialog::getSaveFileName(
                win, "Save image", QString(),
                "PNG image (*.png);;PDF document (*.pdf);;JPEG image (*.jpg)");
            if (path.isEmpty()) return;
            if (path.endsWith(".pdf", Qt::CaseInsensitive))
                customPlot->savePdf(path);
            else if (path.endsWith(".jpg",  Qt::CaseInsensitive) ||
                     path.endsWith(".jpeg", Qt::CaseInsensitive))
                customPlot->saveJpg(path);
            else
                customPlot->savePng(path.endsWith(".png", Qt::CaseInsensitive) ? path : path + ".png");
        });
        toolbar->addWidget(saveBtn);
    }

    win->addToolBar(Qt::TopToolBarArea, toolbar);
}

// Viewer for HDUs that carry no image data (NAXIS=0): shows the header table directly.
static QWidget* makeHeaderOnlyViewer(const HduInfo& hdu, const std::string& fitsPath)
{
    DataFile fm(fitsPath);
    const FITShdu fitshdu = *fm.GetHeaderAtIndex(hdu.index);

    auto* win = new QMainWindow;
    win->setWindowTitle(QString("HDU %1 — %2 (header only)")
        .arg(hdu.index)
        .arg(QString::fromStdString(hdu.extname)));

    win->setCentralWidget(makeFitsHeaderTable(fitshdu, win));
    win->resize(700, 400);
    return win;
}

static void addCatalogueTools(QMainWindow* win, QCustomPlot* plot,
                              const std::shared_ptr<const FITScube>& img, const std::string& fitsPath);

static double scaledPixel(double value, QtFITS::zScale scale)
{
    switch (scale)
    {
        case QtFITS::fLOG:  return std::log10(value);
        case QtFITS::fSQRT: return std::sqrt(value);
        case QtFITS::fQUAD: return value * value;
        default:            return value;
    }
}

// The pixel values the colour range starts from when it isn't the image's whole range: --zscale if
// given, else the percentiles if the image's dynamic is large, as QtFITSviewer::Overlay does.
static bool chosenPixelLimits(const std::shared_ptr<const FITScube>& img, double& lower, double& upper)
{
    if (minZscale == minZscale && maxZscale == maxZscale) {
        lower = minZscale;
        upper = maxZscale;
        return true;
    }
    if (std::abs(img->GetMaximum() / img->GetMinimum()) > 100) {
        lower = img->Get5thpercentil();
        upper = img->Get95thpercentil();
        return true;
    }
    return false;
}

// Sinh stretch of a pixel, from 0 at lower to 1 at upper. sinh(pixel) itself overflows from a pixel
// value of about 710, so the value is first brought to [0, 1]. The factor 3 is the one DS9 uses.
static double sinhStretch(double value, double lower, double upper)
{
    const double position = std::min(1.0, std::max(0.0, (value - lower) / (upper - lower)));
    return std::sinh(3.0 * position) / std::sinh(3.0);
}

// Rewrites the map's cells from the pixels [iMin, iMax[ of the image, through the scale, and names
// the scale on the colour bar. The map is kept rather than redrawn: a new one would hide what is
// overlaid and reset the zoom. A value the scale leaves undefined, such as the log of a negative
// pixel, takes the end of the range it fell off. Returns the range of the cells.
static QCPRange fillScaledImage(QCPColorMap* map, QCPColorScale* colorScale, const std::shared_ptr<const FITScube>& img,
                                size_t iMin, size_t iMax, QtFITS::zScale scale)
{
    const std::valarray<bool> mask = img->GetMask();
    std::vector<double> values(iMax - iMin);

    double sinhLower = 0.0, sinhUpper = 0.0;
    if (scale == QtFITS::fSINH && (!chosenPixelLimits(img, sinhLower, sinhUpper) || !(sinhLower < sinhUpper))) {
        sinhLower = img->GetMinimum();
        sinhUpper = img->GetMaximum();
    }

    double lowest  =  std::numeric_limits<double>::max();
    double highest = -std::numeric_limits<double>::max();
    for (size_t ik = iMin; ik < iMax; ++ik) {
        if (mask[ik])
            continue;
        const double value = (scale == QtFITS::fSINH) ? sinhStretch(img->at(ik), sinhLower, sinhUpper)
                                                      : scaledPixel(img->at(ik), scale);
        values[ik - iMin] = value;
        if (std::isfinite(value)) {
            lowest  = std::min(lowest, value);
            highest = std::max(highest, value);
        }
    }
    if (lowest > highest) {
        lowest  = 0.0;
        highest = 1.0;
    } else if (lowest == highest) {
        highest = lowest + 1.0;
    }

    for (size_t ik = iMin; ik < iMax; ++ik) {
        if (mask[ik])
            continue;
        double value = values[ik - iMin];
        if (!std::isfinite(value))
            value = (value == std::numeric_limits<double>::infinity()) ? highest : lowest;
        const std::vector<size_t> pixel = img->PixelCoordinates(ik);
        map->data()->setCell(static_cast<int>(pixel[0]), static_cast<int>(pixel[1]), value);
    }

    const std::string units = img->HDU().Exists("UNITS") ? img->HDU().GetValueForKey("UNITS") : "ARB. UNIT";
    switch (scale)
    {
        case QtFITS::fLOG:  colorScale->axis()->setLabel(QString::fromStdString("log10(" + units + ")")); break;
        case QtFITS::fSQRT: colorScale->axis()->setLabel(QString::fromStdString(units + "^(1/2)"));       break;
        case QtFITS::fQUAD: colorScale->axis()->setLabel(QString::fromStdString(units + "^2"));           break;
        case QtFITS::fSINH: colorScale->axis()->setLabel(QString::fromStdString("sinh(" + units + "), normalised")); break;
        default:            colorScale->axis()->setLabel(QString::fromStdString(units));                  break;
    }
    return QCPRange(lowest, highest);
}

// The colour range to start from under a scale: the chosen limits if there are any, else the cells'
// whole range. The limits are in pixel values, so they go through the scale too, and fall back on
// the cells' range if the scale leaves them undefined. The sinh stretch has already used them.
static QCPRange initialScaleRange(const std::shared_ptr<const FITScube>& img, QtFITS::zScale scale, const QCPRange& cells)
{
    if (scale == QtFITS::fSINH)
        return QCPRange(0.0, 1.0);

    double lower, upper;
    if (!chosenPixelLimits(img, lower, upper))
        return cells;

    const bool crossesZero = lower < 0.0 && upper > 0.0;
    const double a = scaledPixel(lower, scale);
    const double b = scaledPixel(upper, scale);
    if (!std::isfinite(a) || !std::isfinite(b) || a == b)
        return cells;
    if (scale == QtFITS::fQUAD && crossesZero)
        return QCPRange(0.0, std::max(a, b));
    return QCPRange(std::min(a, b), std::max(a, b));
}

// Adds to the image's toolbar a choice of the display scale, shared through scale with whatever
// redraws the image. pixels gives the pixels of the image displayed, [first, last[.
static void addScaleTool(QMainWindow* win, QCustomPlot* plot, const std::shared_ptr<const FITScube>& img,
                         const std::shared_ptr<QtFITS::zScale>& scale,
                         const std::function<std::pair<size_t, size_t>()>& pixels)
{
    auto* toolbar = win->findChild<QToolBar*>();
    if (!toolbar || !plot->plotLayout()->hasElement(0, 1))
        return;
    auto* colorScale = dynamic_cast<QCPColorScale*>(plot->plotLayout()->element(0, 1));
    if (!colorScale)
        return;

    auto* scaleCombo = new QComboBox(toolbar);
    scaleCombo->setToolTip("Scale the pixel values go through to be displayed");
    scaleCombo->addItem("Linear",    static_cast<int>(QtFITS::fLIN));
    scaleCombo->addItem("Sqrt",      static_cast<int>(QtFITS::fSQRT));
    scaleCombo->addItem("Log10",     static_cast<int>(QtFITS::fLOG));
    scaleCombo->addItem("Quadratic", static_cast<int>(QtFITS::fQUAD));
    scaleCombo->addItem("Sinh",      static_cast<int>(QtFITS::fSINH));

    toolbar->addSeparator();
    toolbar->addWidget(new QLabel("Scale:", toolbar));
    toolbar->addWidget(scaleCombo);

    QObject::connect(scaleCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), win, [=](int) {
        const QList<QCPColorMap*> maps = colorScale->colorMaps();
        if (maps.isEmpty())
            return;
        *scale = static_cast<QtFITS::zScale>(scaleCombo->currentData().toInt());
        const std::pair<size_t, size_t> range = pixels();
        const QCPRange cells = fillScaledImage(maps.last(), colorScale, img, range.first, range.second, *scale);
        colorScale->setDataRange(initialScaleRange(img, *scale, cells));
        plot->replot();
    });
}

// Adds to the image's toolbar a choice of the colour palette the image is drawn with. The colour
// scale holds the palette and passes it on to the image. The first entry is QtDisplay's default.
static void addPaletteTool(QMainWindow* win, QCustomPlot* plot)
{
    auto* toolbar = win->findChild<QToolBar*>();
    if (!toolbar || !plot->plotLayout()->hasElement(0, 1))
        return;
    auto* colorScale = dynamic_cast<QCPColorScale*>(plot->plotLayout()->element(0, 1));
    if (!colorScale)
        return;

    struct Palette { const char* name; QCPColorGradient gradient; };
    const std::vector<Palette> palettes = {
        {"Night",     QCPColorGradient::gpNight},
        {"Grayscale", QCPColorGradient::gpGrayscale},
        {"Viridis",   QtFITS::makeViridisGradient()},
        {"Hot",       QCPColorGradient::gpHot},
        {"Cold",      QCPColorGradient::gpCold},
        {"Thermal",   QCPColorGradient::gpThermal},
        {"Ion",       QCPColorGradient::gpIon},
        {"Candy",     QCPColorGradient::gpCandy},
        {"Geography", QCPColorGradient::gpGeography},
        {"Polar",     QCPColorGradient::gpPolar},
        {"Spectrum",  QCPColorGradient::gpSpectrum},
        {"Jet",       QCPColorGradient::gpJet},
        {"Hues",      QCPColorGradient::gpHues}
    };

    auto* paletteCombo = new QComboBox(toolbar);
    paletteCombo->setToolTip("Colour palette of the image");
    for (const auto& palette : palettes)
        paletteCombo->addItem(palette.name);

    toolbar->addSeparator();
    toolbar->addWidget(new QLabel("Palette:", toolbar));
    toolbar->addWidget(paletteCombo);

    QObject::connect(paletteCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), win, [=](int index) {
        colorScale->setGradient(palettes[static_cast<size_t>(index)].gradient);
        plot->replot();
    });
}

static QWidget* make2DViewer(const HduInfo& hdu,
                              const std::string& fitsPath,
                              const UOption& opt)
{
    // Thread cleanup is now handled by makeViewer before dispatching here.

    DataFile fm(fitsPath);
    std::shared_ptr<FITScube> img = fm.GetImageAtIndex(hdu.index);

    if (opt.ROI.size() > 1) {
        size_t X0 = 0, nX = img->Size(1);
        size_t Y0 = 0, nY = img->Size(2);
        if (opt.ROI.find('x') != opt.ROI.end()) {
            X0 = opt.ROI.find('x')->second.first;
            nX = opt.ROI.find('x')->second.second;
        }
        if (opt.ROI.find('y') != opt.ROI.end()) {
            Y0 = opt.ROI.find('y')->second.first;
            nY = opt.ROI.find('y')->second.second;
        }
        auto tmp = img->Window(X0, Y0, nX, nY);
        std::swap(img, tmp);
    }

    auto c_img = std::const_pointer_cast<const FITScube>(img);
    QCustomPlot* plt = Draw(c_img);
    if (!plt) return nullptr;

    QMainWindow* mainWin = qobject_cast<QMainWindow*>(plt->window());
    if (mainWin) {
        addCommandBar(mainWin, img->HDU(), plt);
        addScaleTool(mainWin, plt, c_img, std::make_shared<QtFITS::zScale>(QtFITS::fLIN),
                     [c_img]() { return std::make_pair(size_t(0), static_cast<size_t>(c_img->Nelements())); });
        addPaletteTool(mainWin, plt);
        addCatalogueTools(mainWin, plt, c_img, fitsPath);
    }
    return mainWin ? mainWin : dynamic_cast<QWidget*>(plt->window());
}

static QWidget* make3DViewer(const HduInfo& hdu,
                              const std::string& fitsPath,
                              const UOption&)
{
    DataFile fm(fitsPath);
    std::shared_ptr<FITScube> cube = fm.GetImageAtIndex(hdu.index);
    auto c_cube = std::const_pointer_cast<const FITScube>(cube);
    int nLayers = hdu.naxis3;
    auto scale = std::make_shared<QtFITS::zScale>(QtFITS::fLIN);

    auto* win = new QMainWindow;
    win->setWindowTitle(QString("HDU %1 - %2 (%3 layers)")
                        .arg(hdu.index)
                        .arg(QString::fromStdString(hdu.extname))
                        .arg(nLayers));

    auto* central = new QWidget(win);
    auto* vlay = new QVBoxLayout(central);
    vlay->setContentsMargins(0, 0, 0, 0);
    vlay->setSpacing(4);

    auto* customPlot = new QCustomPlot(central);
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes);
    customPlot->axisRect()->setupFullAxesBox(true);
    vlay->addWidget(customPlot, 1);

    QtFITSviewer::QtDisplayLayer(customPlot, c_cube, 1);

    QCPAxisRect* image = customPlot->axisRect(0);
    image->setBackgroundScaledMode(Qt::AspectRatioMode::KeepAspectRatioByExpanding);

    QCPLayoutGrid* xproj_layout = new QCPLayoutGrid;
    customPlot->plotLayout()->addElement(1, 0, xproj_layout);
    QCPLayoutGrid* histo_layout = new QCPLayoutGrid;
    customPlot->plotLayout()->addElement(0, 3, histo_layout);
    QCPLayoutGrid* yproj_layout = new QCPLayoutGrid;
    customPlot->plotLayout()->addElement(0, 2, yproj_layout);

    QCPAxisRect* yproj = new QCPAxisRect(customPlot, false);
    QCPAxisRect* xproj = new QCPAxisRect(customPlot, false);
    QCPAxisRect* histo = new QCPAxisRect(customPlot, false);

    image->setMinimumSize(200, 200);
    yproj->setMinimumSize(50, 200);  yproj->setMaximumSize(100, 1200);
    xproj->setMinimumSize(200, 50);  xproj->setMaximumSize(1200, 100);
    histo->setMinimumSize(60, 200);

    yproj_layout->addElement(0, 0, yproj);
    xproj_layout->addElement(0, 0, xproj);
    histo_layout->addElement(0, 0, histo);

    xproj->addAxes(QCPAxis::atBottom | QCPAxis::atLeft);
    xproj->axis(QCPAxis::atBottom)->setRange(image->axis(QCPAxis::atBottom)->range());
    xproj->axis(QCPAxis::atLeft)->ticker()->setTickCount(2);
    xproj->axis(QCPAxis::atBottom)->ticker()->setTickCount(2);
    xproj->axis(QCPAxis::atLeft)->setRange(c_cube->GetMinimum(), c_cube->GetMaximum());

    yproj->addAxes(QCPAxis::atBottom | QCPAxis::atRight);
    yproj->axis(QCPAxis::atRight)->setRange(image->axis(QCPAxis::atLeft)->range());
    yproj->axis(QCPAxis::atRight)->ticker()->setTickCount(2);
    yproj->axis(QCPAxis::atBottom)->ticker()->setTickCount(2);
    yproj->axis(QCPAxis::atBottom)->setRangeReversed(true);
    yproj->axis(QCPAxis::atBottom)->setRange(c_cube->GetMinimum(), c_cube->GetMaximum());

    histo->addAxes(QCPAxis::atBottom | QCPAxis::atRight | QCPAxis::atLeft);
    histo->axis(QCPAxis::atBottom)->setRange(c_cube->GetMinimum()*(1-0.05), c_cube->GetMaximum()*(1+0.05));
    histo->axis(QCPAxis::atBottom)->setScaleType(QCPAxis::stLinear);
    histo->axis(QCPAxis::atBottom)->ticker()->setTickCount(2);
    histo->axis(QCPAxis::atRight)->setRange(-0.05, 1.05);
    histo->axis(QCPAxis::atRight)->setScaleType(QCPAxis::stLinear);
    histo->axis(QCPAxis::atRight)->ticker()->setTickCount(2);
    histo->axis(QCPAxis::atLeft)->setScaleType(QCPAxis::stLogarithmic);
    QSharedPointer<QCPAxisTickerLog> logTicker3d(new QCPAxisTickerLog);
    histo->axis(QCPAxis::atLeft)->setTicker(logTicker3d);

    foreach (QCPAxis *axis, yproj->axes())
        axis->setLayer("axes");
    foreach (QCPAxis *axis, xproj->axes())
        axis->setLayer("axes");
    foreach (QCPAxis *axis, histo->axes())
        axis->setLayer("axes");

    QCPGraph* x_projection = customPlot->addGraph(xproj->axis(QCPAxis::atBottom), xproj->axis(QCPAxis::atLeft));
    x_projection->setLineStyle(QCPGraph::lsStepCenter);
    x_projection->setPen(QPen(QColor("#FFA100"), 1.5));

    QCPGraph* y_projection = customPlot->addGraph(yproj->axis(QCPAxis::atRight), yproj->axis(QCPAxis::atBottom));
    y_projection->setLineStyle(QCPGraph::lsStepCenter);
    y_projection->setPen(QPen(QColor("#FFA100"), 1.5));

    QCPGraph* px_pdf = customPlot->addGraph(histo->axis(QCPAxis::atBottom), histo->axis(QCPAxis::atLeft));
    px_pdf->setLineStyle(QCPGraph::lsStepCenter);
    px_pdf->setPen(QPen(QColor(255, 100, 0)));

    QCPGraph* px_cdf = customPlot->addGraph(histo->axis(QCPAxis::atBottom), histo->axis(QCPAxis::atRight));
    px_cdf->setLineStyle(QCPGraph::lsLine);
    px_cdf->setPen(QPen(QColor(30, 40, 255, 150)));

    if (minZscale == minZscale && maxZscale == maxZscale) {
        if (customPlot->plotLayout()->hasElement(0, 1)) {
            auto* cs = dynamic_cast<QCPColorScale*>(customPlot->plotLayout()->element(0, 1));
            if (cs) cs->setDataRange(QCPRange(minZscale, maxZscale));
        }
    }

    size_t nx = c_cube->Size(1), ny = c_cube->Size(2);
    QtFITSviewer::QtHistogram(customPlot, c_cube, histo_bins,
                               (double)(0), (double)(nx * ny),
                               px_pdf, px_cdf);

    customPlot->connect(image->axis(QCPAxis::atBottom), SIGNAL(rangeChanged(QCPRange)),
                        xproj->axis(QCPAxis::atBottom), SLOT(setRange(QCPRange)));
    customPlot->connect(image->axis(QCPAxis::atLeft), SIGNAL(rangeChanged(QCPRange)),
                        yproj->axis(QCPAxis::atRight), SLOT(setRange(QCPRange)));

    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes | QCP::iSelectOther);
    customPlot->setNoAntialiasingOnDrag(true);
    customPlot->setPlottingHint(QCP::phImmediateRefresh, false);
    image->axis(QCPAxis::atLeft)->setScaleRatio(image->axis(QCPAxis::atBottom), 1.0);
    QObject::connect(customPlot, &QCustomPlot::beforeReplot, [image]() {
        image->axis(QCPAxis::atLeft)->setScaleRatio(image->axis(QCPAxis::atBottom), 1.0);
    });

    auto fMouse = std::bind(mouseHasMove, customPlot, c_cube, x_projection, y_projection);
    QAbstractAnimation::connect(customPlot, &QCustomPlot::mouseMove, fMouse);

    // Slider strip
    auto* strip = new QWidget(central);
    auto* hbox = new QHBoxLayout(strip);
    hbox->setContentsMargins(4, 2, 4, 2);
    auto* sliderLabel = new QLabel("Layer:", strip);
    auto* slider = new QSlider(Qt::Horizontal, strip);
    slider->setRange(0, nLayers - 1);
    slider->setValue(0);
    slider->setTickPosition(QSlider::TicksBelow);
    slider->setTickInterval(std::max(1, nLayers / 20));
    // No base style: QProxyStyle takes ownership of the one it is given, and slider->style() is the
    // application's, which would be deleted with the slider and crash every window left.
    auto* arrowStyle = new ArrowSliderStyle();
    arrowStyle->setParent(slider);
    slider->setStyle(arrowStyle);
    auto* layerLabel = new QLabel(QString("1 / %1").arg(nLayers), strip);
    layerLabel->setMinimumWidth(60);
    hbox->addWidget(sliderLabel);
    hbox->addWidget(slider, 1);
    hbox->addWidget(layerLabel);
    vlay->addWidget(strip, 0);

    win->setCentralWidget(central);
    win->resize(1200, 950);

    QObject::connect(slider, &QSlider::valueChanged,
                     [customPlot, c_cube, px_pdf, px_cdf, layerLabel, nLayers, nx, ny, scale](int layer0) {
                         int layer1 = layer0 + 1;

                         // Preserve current color scale range, palette and axis limits across layer
                         // changes: QtDisplayLayer resets the palette to its default
                         QCPRange savedRange, savedXRange, savedYRange;
                         QCPColorGradient savedGradient;
                         bool hasRange = false;
                         if (customPlot->plotLayout()->hasElement(0, 1)) {
                             auto* cs = dynamic_cast<QCPColorScale*>(
                                 customPlot->plotLayout()->element(0, 1));
                             if (cs) { savedRange = cs->dataRange(); savedGradient = cs->gradient(); hasRange = true; }
                         }
                         QCPAxisRect* imageRect = customPlot->axisRect(0);
                         savedXRange = imageRect->axis(QCPAxis::atBottom)->range();
                         savedYRange = imageRect->axis(QCPAxis::atLeft)->range();

                         QtFITSviewer::QtDisplayLayer(customPlot, c_cube, static_cast<size_t>(layer1));

                         // Apply the scale chosen and restore color scale range and palette
                         if (hasRange) {
                             if (customPlot->plotLayout()->hasElement(0, 1)) {
                                 auto* cs = dynamic_cast<QCPColorScale*>(
                                     customPlot->plotLayout()->element(0, 1));
                                 // QtDisplayLayer adds a map per call: the older ones are hidden
                                 // but would still be drawn at every zoom and pan.
                                 if (cs) {
                                     const QList<QCPColorMap*> maps = cs->colorMaps();
                                     for (int i = 0; i + 1 < maps.size(); ++i)
                                         customPlot->removePlottable(maps[i]);
                                 }
                                 if (cs && *scale != QtFITS::fLIN && !cs->colorMaps().isEmpty())
                                     fillScaledImage(cs->colorMaps().last(), cs, c_cube,
                                                     static_cast<size_t>(layer0) * nx * ny,
                                                     static_cast<size_t>(layer1) * nx * ny, *scale);
                                 if (cs) {
                                     cs->setDataRange(savedRange);
                                     cs->setGradient(savedGradient);
                                 }
                             }
                         }
                         // Restore axis limits
                         imageRect->axis(QCPAxis::atBottom)->setRange(savedXRange);
                         imageRect->axis(QCPAxis::atLeft)->setRange(savedYRange);

                         QtFITSviewer::QtHistogram(customPlot, c_cube, histo_bins,
                                                    (double)(layer0) * (double)(nx * ny),
                                                    (double)(layer1) * (double)(nx * ny),
                                                    px_pdf, px_cdf);
                         layerLabel->setText(QString("%1 / %2").arg(layer1).arg(nLayers));
                         customPlot->replot(QCustomPlot::rpQueuedReplot);
                     });

    addCommandBar(win, cube->HDU(), customPlot);
    addScaleTool(win, customPlot, c_cube, scale, [slider, nx, ny]() {
        const size_t layer = static_cast<size_t>(slider->value());
        return std::make_pair(layer * nx * ny, (layer + 1) * nx * ny);
    });
    addPaletteTool(win, customPlot);
    return win;
}

// ---------------------------------------------------------------------------
// FITSTable viewer
// ---------------------------------------------------------------------------

// Returns true if the column can be used as Y axis (ordinate).
// Complex types (tcplx, tdbcplx) are allowed — the plotting backend handles them.
static bool isPlottableAsY(const std::unique_ptr<FITSform>& col)
{
    const DSL::dtype t = col->getType();
    return t != DSL::dtype::tsbyte
        && t != DSL::dtype::tstring
        && t != DSL::dtype::tlogical
        && t != DSL::dtype::tbit
        && t != DSL::dtype::tbyte
        && t != DSL::dtype::tnone;
}

// Returns true if the column can be used as X axis (abscissa).
// Complex types are excluded: no natural ordering as a real-valued axis.
static bool isPlottableAsX(const std::unique_ptr<FITSform>& col)
{
    return isPlottableAsY(col)
        && col->getType() != DSL::dtype::tcplx
        && col->getType() != DSL::dtype::tdbcplx;
}

// Rows of a table sharing one category of the column chosen to split a plot.
struct TableCategory { QString label; RowSet rows; QColor color{}; };

// Keeps a continuous column from giving one graph per row. Colors can't identify this many
// categories by eye: clicking one in the plot's category list highlights it.
static const size_t kMaxTableCategories = 512;

// Above this, evenly spaced hues stop being distinguishable.
static const size_t kMaxHueCategories = 32;

// Above this, a color ramp no longer separates neighbouring numeric categories.
static const size_t kMaxRampCategories = 360;

// Returns true if the column can split a plot: any scalar column except complex, or a multi-bit column.
static bool isCategoryColumn(const FITSform& col)
{
    const std::type_index payload = col.payloadType();
    return payload == typeid(int8_t)   || payload == typeid(uint8_t)
        || payload == typeid(int16_t)  || payload == typeid(uint16_t)
        || payload == typeid(int32_t)  || payload == typeid(uint32_t)
        || payload == typeid(int64_t)  || payload == typeid(uint64_t)
        || payload == typeid(float)    || payload == typeid(double)
        || payload == typeid(bool)     || payload == typeid(std::string)
        || (payload == typeid(FITSform::boolVector) && col.getType() == DSL::dtype::tbit);
}

static bool isRealColumn(const FITSform& col)
{
    return col.getType() == DSL::dtype::tfloat || col.getType() == DSL::dtype::tdouble;
}

// Splits rows by the key keyOf() gives them; categories come out in ascending key order.
template<typename Key, typename KeyFn, typename LabelFn>
static std::vector<TableCategory> groupRows(const RowSet& rows, KeyFn keyOf, LabelFn labelOf)
{
    std::map<Key, std::vector<size_t>> groups;
    for (size_t row : rows.indices()) {
        groups[keyOf(row)].push_back(row);
        if (groups.size() > kMaxTableCategories)
            throw std::runtime_error("More than " + std::to_string(kMaxTableCategories)
                + " categories in this column. For a float or double column, try a larger epsilon.");
    }

    std::vector<TableCategory> categories;
    categories.reserve(groups.size());
    for (auto& group : groups)
        categories.push_back({labelOf(group.first), RowSet(std::move(group.second))});
    return categories;
}

// One category per distinct value (integer, byte, logical, single bit, string).
template<typename T>
static std::vector<TableCategory> categoriesOfValues(const FITSform& col, const RowSet& rows)
{
    const auto& values = col.values<T>();
    return groupRows<T>(rows,
        [&values](size_t row) -> T { return values[row]; },
        [](const T& value) -> QString {
            if constexpr (std::is_same_v<T, std::string>) return QString::fromStdString(value);
            else if constexpr (std::is_same_v<T, bool>)   return value ? "T" : "F";
            else                                          return QString::number(value);
        });
}

// One category per contiguous interval of width 2*epsilon, starting with [v0-epsilon, v0+epsilon]
// around the column's smallest finite value v0: every value falls in exactly one, val+epsilon
// still belongs to the category of val, and anything above it to the next one.
// v0 comes from the whole column so that an interval doesn't move with the rows selected.
// With epsilon = 0, one category per distinct value. NaN and infinities get one category each.
template<typename T>
static std::vector<TableCategory> categoriesOfReals(const FITSform& col, const RowSet& rows, double epsilon)
{
    const auto& values = col.values<T>();

    double v0 = std::numeric_limits<double>::infinity();
    for (const T value : values)
        if (std::isfinite(value))
            v0 = std::min(v0, static_cast<double>(value));

    // {rank, position}: rank orders -Inf, finite values, +Inf, NaN.
    using Key = std::pair<int, double>;
    enum { kMinusInf, kFinite, kPlusInf, kNaN };

    const int digits = std::numeric_limits<T>::digits10;
    const QString plusMinus = QString(QChar(0x00B1));

    return groupRows<Key>(rows,
        [&values, v0, epsilon](size_t row) -> Key {
            const double value = static_cast<double>(values[row]);
            if (std::isnan(value)) return {kNaN, 0.0};
            if (std::isinf(value)) return {value < 0.0 ? kMinusInf : kPlusInf, 0.0};
            if (epsilon == 0.0)    return {kFinite, value};
            return {kFinite, std::max(0.0, std::ceil((value - v0 - epsilon) / (2.0 * epsilon)))};
        },
        [v0, epsilon, digits, plusMinus](const Key& key) -> QString {
            switch (key.first) {
                case kNaN:      return "NaN";
                case kMinusInf: return "-Inf";
                case kPlusInf:  return "+Inf";
                default:        break;
            }
            if (epsilon == 0.0)
                return QString::number(key.second, 'g', digits);
            return QString("%1 %2 %3")
                .arg(QString::number(v0 + 2.0 * epsilon * key.second, 'g', digits))
                .arg(plusMinus)
                .arg(QString::number(epsilon, 'g', digits));
        });
}

// One category per bit pattern, e.g. 01001010.
static std::vector<TableCategory> categoriesOfBits(const FITSform& col, const RowSet& rows)
{
    const auto& values = col.values<FITSform::boolVector>();
    return groupRows<std::string>(rows,
        [&values](size_t row) -> std::string {
            std::string pattern;
            pattern.reserve(values[row].size());
            for (const bool bit : values[row])
                pattern += bit ? '1' : '0';
            return pattern;
        },
        [](const std::string& pattern) -> QString { return QString::fromStdString(pattern); });
}

// Splits rows by the column's values, in ascending order. Throws above kMaxTableCategories.
static std::vector<TableCategory> groupByColumn(const FITSform& col, const RowSet& rows, double epsilon)
{
    const std::type_index payload = col.payloadType();

    if (payload == typeid(double))               return categoriesOfReals<double>(col, rows, epsilon);
    if (payload == typeid(float))                return categoriesOfReals<float>(col, rows, epsilon);
    if (payload == typeid(int8_t))               return categoriesOfValues<int8_t>(col, rows);
    if (payload == typeid(uint8_t))              return categoriesOfValues<uint8_t>(col, rows);
    if (payload == typeid(int16_t))              return categoriesOfValues<int16_t>(col, rows);
    if (payload == typeid(uint16_t))             return categoriesOfValues<uint16_t>(col, rows);
    if (payload == typeid(int32_t))              return categoriesOfValues<int32_t>(col, rows);
    if (payload == typeid(uint32_t))             return categoriesOfValues<uint32_t>(col, rows);
    if (payload == typeid(int64_t))              return categoriesOfValues<int64_t>(col, rows);
    if (payload == typeid(uint64_t))             return categoriesOfValues<uint64_t>(col, rows);
    if (payload == typeid(bool))                 return categoriesOfValues<bool>(col, rows);
    if (payload == typeid(std::string))          return categoriesOfValues<std::string>(col, rows);
    if (payload == typeid(FITSform::boolVector)) return categoriesOfBits(col, rows);

    throw std::invalid_argument("Column '" + col.getName() + "' can't be used as a category.");
}

// Numeric categories get evenly spaced hues while these can still tell them apart, then follow
// Viridis in ascending order, so that the colors show the trend. Others, and numeric ones above
// kMaxRampCategories, step by the golden angle, so that neighbouring categories get clearly different hues.
static QColor categoryColor(size_t index, size_t count, bool numeric)
{
    if (numeric && count <= kMaxHueCategories)
        return QColor::fromHsvF(static_cast<float>(index) / static_cast<float>(count), 0.90f, 0.78f);

    if (numeric && count <= kMaxRampCategories) {
        static QCPColorGradient ramp = [] {
            QCPColorGradient gradient = QtFITS::makeViridisGradient();
            gradient.setLevelCount(1024); // the default 350 levels would merge categories
            return gradient;
        }();
        // Stops short of Viridis' last yellows, too faint on the plot's white background.
        const double position = count > 1 ? 0.9 * static_cast<double>(index) / static_cast<double>(count - 1) : 0.0;
        return QColor::fromRgb(ramp.color(position, QCPRange(0.0, 1.0)));
    }

    const double goldenAngleInTurns = 0.381966011250105;
    const double hue = std::fmod(static_cast<double>(index) * goldenAngleInTurns, 1.0);
    return QColor::fromHsvF(static_cast<float>(hue), 0.90f, 0.78f);
}

// Splits rows by the column's values and gives each category its color.
static std::vector<TableCategory> buildCategories(const FITSform& col, const RowSet& rows, double epsilon)
{
    std::vector<TableCategory> categories = groupByColumn(col, rows, epsilon);

    const std::type_index payload = col.payloadType();
    const bool numeric = payload != typeid(std::string) && payload != typeid(bool)
                      && payload != typeid(FITSform::boolVector);
    for (size_t i = 0; i < categories.size(); ++i)
        categories[i].color = categoryColor(i, categories.size(), numeric);
    return categories;
}

// ---------------------------------------------------------------------------
// Copying from the table views
// ---------------------------------------------------------------------------

// Copies the selected cells as text: tab between columns, new line between rows. A model
// can give the full text through Qt::UserRole when what it displays is shortened.
static void copySelectedCells(QAbstractItemView* view)
{
    // Every cell becomes text at once here, unlike on screen: "select all" on a huge column stops here.
    static const qint64 kMaxCopiedCells = 1000000;

    qint64 cells = 0;
    for (const QItemSelectionRange& range : view->selectionModel()->selection())
        cells += static_cast<qint64>(range.width()) * range.height();

    if (cells > kMaxCopiedCells) {
        QMessageBox::warning(view, "Copy", QString("%1 cells are selected: too many to copy, the limit is %2.")
            .arg(cells).arg(kMaxCopiedCells));
        return;
    }

    QModelIndexList indexes = view->selectionModel()->selectedIndexes();
    if (indexes.isEmpty())
        return;
    std::sort(indexes.begin(), indexes.end()); // by row, then by column

    QString text;
    for (int i = 0; i < indexes.size(); ++i) {
        if (i > 0)
            text += QLatin1Char(indexes[i].row() == indexes[i - 1].row() ? '\t' : '\n');
        const QVariant full = indexes[i].data(Qt::UserRole);
        text += (full.isValid() ? full : indexes[i].data(Qt::DisplayRole)).toString();
    }
    QApplication::clipboard()->setText(text);
}

// The usual copy key on a table view: on its own, Qt's only copies the current cell.
class CellCopyKeyFilter : public QObject
{
public:
    explicit CellCopyKeyFilter(QAbstractItemView* view) : QObject(view), m_view(view) {}

    bool eventFilter(QObject*, QEvent* e) override
    {
        if (e->type() == QEvent::KeyPress && static_cast<QKeyEvent*>(e)->matches(QKeySequence::Copy)) {
            copySelectedCells(m_view);
            return true;
        }
        return false;
    }

private:
    QAbstractItemView* m_view;
};

// Copy key and right-click menu on the cells of a view.
static void enableCellCopy(QAbstractItemView* view)
{
    view->installEventFilter(new CellCopyKeyFilter(view));

    view->setContextMenuPolicy(Qt::CustomContextMenu);
    QObject::connect(view, &QWidget::customContextMenuRequested, view, [view](const QPoint& position) {
        QMenu menu(view);
        QObject::connect(menu.addAction("Copy"), &QAction::triggered, view, [view]() { copySelectedCells(view); });
        menu.exec(view->viewport()->mapToGlobal(position));
    });
}

// Right-click menu on the column headers of a view, copying the column's name.
static void enableColumnNameCopy(QTableView* view)
{
    QHeaderView* header = view->horizontalHeader();
    header->setContextMenuPolicy(Qt::CustomContextMenu);
    QObject::connect(header, &QWidget::customContextMenuRequested, view, [view, header](const QPoint& position) {
        const int section = header->logicalIndexAt(position);
        if (section < 0)
            return;

        // Qt::UserRole gives the bare name, when the header also shows a unit or an element.
        const QVariant name = view->model()->headerData(section, Qt::Horizontal, Qt::UserRole);
        const QString text = (name.isValid() ? name : view->model()->headerData(section, Qt::Horizontal)).toString();

        QMenu menu(view);
        QObject::connect(menu.addAction("Copy column name"), &QAction::triggered, view,
                         [text]() { QApplication::clipboard()->setText(text); });
        menu.exec(header->mapToGlobal(position));
    });
}

// ---------------------------------------------------------------------------
// Column window: values of one table column, with its properties and statistics
// ---------------------------------------------------------------------------

template<typename T>
static QString cellText(const T& value)
{
    if constexpr (std::is_same_v<T, std::string>)       return QString::fromStdString(value);
    else if constexpr (std::is_same_v<T, bool>)         return value ? "T" : "F";
    else if constexpr (std::is_floating_point_v<T>)     return QString::number(value, 'g', std::numeric_limits<T>::digits10);
    else                                                return QString::number(value);
}

template<typename T>
static QString cellText(const std::pair<T, T>& value)
{
    return QString("%1 %2 %3i")
        .arg(cellText(value.first))
        .arg(QLatin1Char(value.second < T(0) ? '-' : '+'))
        .arg(cellText(static_cast<T>(std::abs(value.second))));
}

// Formats one cell on request, reading straight from the column's own storage.
struct ColumnCellReader
{
    std::function<QString(size_t row, size_t element)> read;
    int  elements     = 1;
    bool rightAligned = true;
};

template<typename T>
static ColumnCellReader scalarCellReader(const FITSform& col)
{
    const auto* values = &col.values<T>();
    ColumnCellReader reader;
    reader.read = [values](size_t row, size_t) -> QString {
        const T value = (*values)[row];
        return cellText(value);
    };
    reader.rightAligned = !std::is_same_v<T, std::string>;
    return reader;
}

template<typename T>
static ColumnCellReader vectorCellReader(const FITSform& col)
{
    const auto* values = &col.values<std::vector<T>>();
    ColumnCellReader reader;
    reader.read = [values](size_t row, size_t element) -> QString {
        const auto& cell = (*values)[row];
        if (element >= cell.size())
            return QString();
        const T value = cell[element];
        return cellText(value);
    };
    reader.elements = static_cast<int>(std::clamp<int64_t>(col.getNelem(), 1, std::numeric_limits<int>::max()));
    reader.rightAligned = !std::is_same_v<T, std::string>;
    return reader;
}

// A multi-bit column reads best as one pattern, e.g. 01001010.
static ColumnCellReader bitPatternCellReader(const FITSform& col)
{
    const auto* values = &col.values<FITSform::boolVector>();
    ColumnCellReader reader;
    reader.read = [values](size_t row, size_t) -> QString {
        QString pattern;
        pattern.reserve(static_cast<int>((*values)[row].size()));
        for (const bool bit : (*values)[row])
            pattern += QLatin1Char(bit ? '1' : '0');
        return pattern;
    };
    return reader;
}

template<typename... Ts>
static ColumnCellReader cellReaderOfPayload(const FITSform& col)
{
    const std::type_index payload = col.payloadType();
    ColumnCellReader reader;
    ((payload == typeid(Ts)              ? void(reader = scalarCellReader<Ts>(col))
    : payload == typeid(std::vector<Ts>) ? void(reader = vectorCellReader<Ts>(col))
    : void()), ...);
    return reader;
}

static ColumnCellReader makeColumnCellReader(const FITSform& col)
{
    if (col.getType() == DSL::dtype::tbit && col.payloadType() == typeid(FITSform::boolVector))
        return bitPatternCellReader(col);

    return cellReaderOfPayload<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t, int64_t, uint64_t,
                               float, double, bool, std::string,
                               FITSform::complex, FITSform::dblcomplex>(col);
}

// QTableView only asks for the cells on screen, and only while they are on screen. This model
// answers from the column's storage and keeps nothing: no copy of the values and no item per
// cell, so scrolling costs the same whatever the size of the column.
class FITSColumnModel : public QAbstractTableModel
{
public:
    FITSColumnModel(std::shared_ptr<const FITStable> table, const FITSform& col, QObject* parent)
        : QAbstractTableModel(parent)
        , m_table(std::move(table))
        , m_reader(makeColumnCellReader(col))
        , m_rows(static_cast<int>(std::min<size_t>(col.size(), std::numeric_limits<int>::max())))
        , m_name(QString::fromStdString(col.getName()))
    {
    }

    int rowCount(const QModelIndex& parent = QModelIndex()) const override
    {
        return parent.isValid() ? 0 : m_rows;
    }

    int columnCount(const QModelIndex& parent = QModelIndex()) const override
    {
        return parent.isValid() ? 0 : m_reader.elements;
    }

    QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override
    {
        if (!index.isValid() || !m_reader.read)
            return QVariant();
        if (role == Qt::DisplayRole)
            return m_reader.read(static_cast<size_t>(index.row()), static_cast<size_t>(index.column()));
        if (role == Qt::TextAlignmentRole)
            return static_cast<int>((m_reader.rightAligned ? Qt::AlignRight : Qt::AlignLeft) | Qt::AlignVCenter);
        return QVariant();
    }

    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override
    {
        if (role == Qt::UserRole && orientation == Qt::Horizontal)
            return m_name;
        if (role != Qt::DisplayRole)
            return QVariant();
        if (orientation == Qt::Vertical)
            return section + 1; // FITS rows count from 1, as the plots' "Row index" axis does
        return m_reader.elements > 1 ? QString("[%1]").arg(section) : m_name;
    }

private:
    std::shared_ptr<const FITStable> m_table; // keeps the storage m_reader points into alive
    ColumnCellReader m_reader;
    int              m_rows;
    QString          m_name;
};

struct ColumnStats
{
    size_t  valid     = 0;
    size_t  nonFinite = 0;
    double  sum       = 0.0;
    double  min       = 0.0;
    double  max       = 0.0;
    double  mean      = 0.0;
    double  stddev    = 0.0; // sample standard deviation
    double  median    = 0.0;
    QString error;
};

// Logical and bit values count as 0 and 1, so their sum is the number of true values.
template<typename T>
static void appendColumnValues(const FITSform& col, std::vector<double>& out)
{
    const std::type_index payload = col.payloadType();
    if (payload == typeid(T)) {
        for (const T value : col.values<T>())
            out.push_back(static_cast<double>(value));
    } else if (payload == typeid(std::vector<T>)) {
        for (const auto& cell : col.values<std::vector<T>>())
            for (const T value : cell)
                out.push_back(static_cast<double>(value));
    }
}

template<typename T>
static void appendColumnModuli(const FITSform& col, std::vector<double>& out)
{
    const std::type_index payload = col.payloadType();
    if (payload == typeid(std::pair<T, T>)) {
        for (const auto& value : col.values<std::pair<T, T>>())
            out.push_back(std::hypot(static_cast<double>(value.first), static_cast<double>(value.second)));
    } else if (payload == typeid(std::vector<std::pair<T, T>>)) {
        for (const auto& cell : col.values<std::vector<std::pair<T, T>>>())
            for (const auto& value : cell)
                out.push_back(std::hypot(static_cast<double>(value.first), static_cast<double>(value.second)));
    }
}

// The median needs a copy of the values: it lives only for the duration of this call.
static ColumnStats computeColumnStats(const FITSform& col)
{
    ColumnStats stats;
    try {
        std::vector<double> values;
        values.reserve(col.size() * static_cast<size_t>(std::max<int64_t>(col.getNelem(), 1)));

        appendColumnValues<int8_t>(col, values);
        appendColumnValues<uint8_t>(col, values);
        appendColumnValues<int16_t>(col, values);
        appendColumnValues<uint16_t>(col, values);
        appendColumnValues<int32_t>(col, values);
        appendColumnValues<uint32_t>(col, values);
        appendColumnValues<int64_t>(col, values);
        appendColumnValues<uint64_t>(col, values);
        appendColumnValues<float>(col, values);
        appendColumnValues<double>(col, values);
        appendColumnValues<bool>(col, values);
        appendColumnModuli<float>(col, values);
        appendColumnModuli<double>(col, values);

        const size_t total = values.size();
        values.erase(std::remove_if(values.begin(), values.end(),
                                    [](double value) { return !std::isfinite(value); }),
                     values.end());
        stats.valid     = values.size();
        stats.nonFinite = total - values.size();
        if (values.empty())
            return stats;

        const auto range = std::minmax_element(values.begin(), values.end());
        stats.min = *range.first;
        stats.max = *range.second;

        long double sum = 0.0L;
        for (const double value : values)
            sum += value;
        stats.sum  = static_cast<double>(sum);
        stats.mean = static_cast<double>(sum / static_cast<long double>(values.size()));

        if (values.size() > 1) {
            long double squares = 0.0L;
            for (const double value : values)
                squares += (static_cast<long double>(value) - stats.mean) * (static_cast<long double>(value) - stats.mean);
            stats.stddev = static_cast<double>(std::sqrt(squares / static_cast<long double>(values.size() - 1)));
        }

        const auto middle = values.begin() + static_cast<std::ptrdiff_t>(values.size() / 2);
        std::nth_element(values.begin(), middle, values.end());
        stats.median = *middle;
        if (values.size() % 2 == 0)
            stats.median = 0.5 * (stats.median + *std::max_element(values.begin(), middle));
    } catch (const std::exception& e) {
        stats.error = QString::fromStdString(e.what());
    }
    return stats;
}

// Adds "name: value" fields to a grid, three per line; returns the label holding the value.
static QLabel* addColumnField(QGridLayout* grid, const QString& name, const QString& value)
{
    const int position = grid->count() / 2;
    auto* valueLabel = new QLabel(value);
    valueLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);
    grid->addWidget(new QLabel(name + ":"), position / 3, 2 * (position % 3));
    grid->addWidget(valueLabel, position / 3, 2 * (position % 3) + 1);
    return valueLabel;
}

// Opens the window of a table column, or brings it forward if it is already open.
// The window shows at once: the statistics are computed in the background and fill in when done.
static void showColumnWindow(QWidget* parent, const std::shared_ptr<const FITStable>& tab,
                             const std::string& name, const std::string& tform)
{
    const QString objectName = "column:" + QString::fromStdString(name);
    if (auto* open = parent->findChild<QWidget*>(objectName, Qt::FindDirectChildrenOnly)) {
        open->raise();
        open->activateWindow();
        return;
    }

    const FITSform& col = *tab->getColumn(name);

    auto* window = new QWidget(parent, Qt::Window);
    window->setObjectName(objectName);
    window->setAttribute(Qt::WA_DeleteOnClose);
    window->setWindowTitle(QString("Column %1 - %2")
        .arg(QString::fromStdString(name), QString::fromStdString(tab->GetName())));
    window->resize(620, 680);
    auto* layout = new QVBoxLayout(window);

    // --- Properties ---
    auto* propertiesBox  = new QGroupBox("Column");
    auto* propertiesGrid = new QGridLayout(propertiesBox);
    addColumnField(propertiesGrid, "Name", QString::fromStdString(name));
    addColumnField(propertiesGrid, "Type", QString::fromStdString(FITSform::getDataType(col.getType())));
    addColumnField(propertiesGrid, "TFORM", QString::fromStdString(tform));
    addColumnField(propertiesGrid, "Unit", col.getUnit().empty() ? QString("-") : QString::fromStdString(col.getUnit()));
    addColumnField(propertiesGrid, "Rows", QString::number(col.size()));
    addColumnField(propertiesGrid, "Elements per row", QString::number(std::max<int64_t>(col.getNelem(), 1)));
    layout->addWidget(propertiesBox);

    // --- Statistics ---
    const std::type_index payload = col.payloadType();
    const bool isText    = payload == typeid(std::string) || payload == typeid(FITSform::stringVector);
    const bool isLogical = col.getType() == DSL::dtype::tlogical || col.getType() == DSL::dtype::tbit;
    const bool isComplex = col.getType() == DSL::dtype::tcplx || col.getType() == DSL::dtype::tdbcplx;
    const bool isVector  = col.getNelem() > 1;

    auto* statsBox  = new QGroupBox(isComplex ? "Statistics (modulus)"
                                  : isVector  ? "Statistics (all elements)"
                                              : "Statistics");
    auto* statsGrid = new QGridLayout(statsBox);
    layout->addWidget(statsBox);

    if (isText) {
        statsGrid->addWidget(new QLabel("No statistics for a string column."), 0, 0);
    } else {
        const QString pending = "computing...";
        std::vector<QLabel*> fields;
        if (isLogical) {
            fields.push_back(addColumnField(statsGrid, "True", pending));
            fields.push_back(addColumnField(statsGrid, "False", pending));
        } else {
            for (const char* field : {"Values", "Non-finite", "Min", "Max", "Mean", "Std dev", "Median"})
                fields.push_back(addColumnField(statsGrid, field, pending));
        }

        // The window can close before the result comes: tab keeps the column alive, and the
        // watcher, owned by the window, goes with it so that nothing updates a deleted label.
        auto* watcher = new QFutureWatcher<ColumnStats>(window);
        QObject::connect(watcher, &QFutureWatcher<ColumnStats>::finished, window, [=]() {
            const ColumnStats stats = watcher->result();

            const auto number = [&stats](double value) {
                return stats.valid == 0 ? QString("-") : QString::number(value, 'g', 10);
            };
            const std::vector<QString> texts = isLogical
                ? std::vector<QString>{QString::number(static_cast<qulonglong>(stats.sum)),
                                       QString::number(static_cast<qulonglong>(stats.valid - static_cast<size_t>(stats.sum)))}
                : std::vector<QString>{QString::number(stats.valid), QString::number(stats.nonFinite),
                                       number(stats.min), number(stats.max), number(stats.mean),
                                       number(stats.stddev), number(stats.median)};

            for (size_t i = 0; i < fields.size(); ++i)
                fields[i]->setText(stats.error.isEmpty() ? texts[i] : QString("failed: ") + stats.error);
        });
        watcher->setFuture(QtConcurrent::run([tab, &col]() { return computeColumnStats(col); }));
    }

    // --- Values ---
    // Fixed row heights and no fit to contents: the view then never measures or visits
    // rows that are off screen, whatever the size of the column.
    auto* view = new QTableView;
    view->verticalHeader()->setSectionResizeMode(QHeaderView::Fixed);
    view->verticalHeader()->setDefaultSectionSize(view->fontMetrics().height() + 6);
    view->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
    view->horizontalHeader()->setDefaultSectionSize(110);
    view->setWordWrap(false);
    view->setAlternatingRowColors(true);
    view->setEditTriggers(QAbstractItemView::NoEditTriggers);

    auto* model = new FITSColumnModel(tab, col, window);
    view->horizontalHeader()->setStretchLastSection(model->columnCount() == 1);
    view->setModel(model);
    enableCellCopy(view);
    enableColumnNameCopy(view);
    layout->addWidget(view, 1);

    window->show();
}

// ---------------------------------------------------------------------------
// Selected rows window: chosen rows of a table, with all its columns
// ---------------------------------------------------------------------------

// Same principle as FITSColumnModel: cells are formatted on request and nothing is kept.
class SelectedRowsModel : public QAbstractTableModel
{
public:
    SelectedRowsModel(std::shared_ptr<const FITStable> table, const RowSet& rows, QObject* parent)
        : QAbstractTableModel(parent)
        , m_table(std::move(table))
        , m_rows(rows.indices())
    {
        for (const auto& entry : m_table->listColumns()) {
            m_readers.push_back(makeColumnCellReader(*m_table->getColumn(entry[0])));
            m_names.push_back(QString::fromStdString(entry[0]));
            m_labels.push_back(m_names.back()
                + (entry[2].empty() ? QString() : " [" + QString::fromStdString(entry[2]) + "]"));
        }
    }

    int rowCount(const QModelIndex& parent = QModelIndex()) const override
    {
        return parent.isValid() ? 0 : static_cast<int>(std::min<size_t>(m_rows.size(), std::numeric_limits<int>::max()));
    }

    int columnCount(const QModelIndex& parent = QModelIndex()) const override
    {
        return parent.isValid() ? 0 : static_cast<int>(m_readers.size());
    }

    QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override
    {
        if (!index.isValid())
            return QVariant();

        const ColumnCellReader& reader = m_readers[static_cast<size_t>(index.column())];
        if (!reader.read)
            return QVariant();

        if (role == Qt::TextAlignmentRole)
            return static_cast<int>((reader.rightAligned && reader.elements == 1 ? Qt::AlignRight : Qt::AlignLeft) | Qt::AlignVCenter);
        if (role != Qt::DisplayRole && role != Qt::UserRole)
            return QVariant();

        const size_t row = m_rows[static_cast<size_t>(index.row())];
        if (reader.elements == 1)
            return reader.read(row, 0);

        // Only the first elements on screen, as a long vector would make the row
        // unreadable; Qt::UserRole, which copying uses, gives them all.
        static const int kShownElements = 6;
        const int nShown = (role == Qt::UserRole) ? reader.elements : std::min(reader.elements, kShownElements);
        QStringList shown;
        for (int e = 0; e < nShown; ++e)
            shown << reader.read(row, static_cast<size_t>(e));
        if (nShown < reader.elements)
            shown << "...";
        return shown.join(", ");
    }

    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override
    {
        if (role == Qt::UserRole && orientation == Qt::Horizontal)
            return m_names[static_cast<size_t>(section)];
        if (role != Qt::DisplayRole)
            return QVariant();
        if (orientation == Qt::Vertical)
            return QVariant::fromValue(static_cast<qulonglong>(m_rows[static_cast<size_t>(section)] + 1));
        return m_labels[static_cast<size_t>(section)];
    }

private:
    std::shared_ptr<const FITStable> m_table; // keeps the storage the readers point into alive
    std::vector<size_t>              m_rows;
    std::vector<ColumnCellReader>    m_readers;
    std::vector<QString>             m_names;  // bare column names
    std::vector<QString>             m_labels; // with the unit, for the headers
};

// Shows rows in the parent's Selected rows window, opened if needed and replaced otherwise.
static void showSelectedRows(QWidget* parent, const std::shared_ptr<const FITStable>& table,
                             const RowSet& rows, const QString& description)
{
    static const char* const kObjectName = "selectedRows";

    QWidget*    window = parent->findChild<QWidget*>(kObjectName, Qt::FindDirectChildrenOnly);
    QTableView* view   = nullptr;

    if (window) {
        view = window->findChild<QTableView*>();
    } else {
        window = new QWidget(parent, Qt::Window);
        window->setObjectName(kObjectName);
        window->setAttribute(Qt::WA_DeleteOnClose);
        window->resize(900, 400);

        view = new QTableView;
        view->verticalHeader()->setSectionResizeMode(QHeaderView::Fixed);
        view->verticalHeader()->setDefaultSectionSize(view->fontMetrics().height() + 6);
        view->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
        view->horizontalHeader()->setDefaultSectionSize(110);
        view->setWordWrap(false);
        view->setAlternatingRowColors(true);
        view->setEditTriggers(QAbstractItemView::NoEditTriggers);
        enableCellCopy(view);
        enableColumnNameCopy(view);

        auto* layout = new QVBoxLayout(window);
        layout->addWidget(view);
    }

    window->setWindowTitle(QString("Selected rows: %1 (%2)").arg(rows.size()).arg(description));

    QAbstractItemModel* previous = view->model();
    view->setModel(new SelectedRowsModel(table, rows, view));
    if (previous)
        previous->deleteLater();

    window->show();
    window->raise();
}

// ---------------------------------------------------------------------------
// Rows behind a plot: scans the plotted columns again, as QCPGraph sorts its
// points by x and so loses the row each one came from
// ---------------------------------------------------------------------------

// All the values of a numeric column as doubles, flattened row after row.
struct ColumnValues
{
    std::vector<double> re;
    std::vector<double> im; // only filled for a complex column
    size_t              elements = 1;
    bool                isComplex = false;
};

template<typename T>
static void appendColumnParts(const FITSform& col, ColumnValues& out)
{
    const std::type_index payload = col.payloadType();
    if (payload == typeid(std::pair<T, T>)) {
        for (const auto& value : col.values<std::pair<T, T>>()) {
            out.re.push_back(static_cast<double>(value.first));
            out.im.push_back(static_cast<double>(value.second));
        }
        out.isComplex = true;
    } else if (payload == typeid(std::vector<std::pair<T, T>>)) {
        for (const auto& cell : col.values<std::vector<std::pair<T, T>>>())
            for (const auto& value : cell) {
                out.re.push_back(static_cast<double>(value.first));
                out.im.push_back(static_cast<double>(value.second));
            }
        out.isComplex = true;
    }
}

static ColumnValues readColumnValues(const FITSform& col)
{
    ColumnValues values;
    values.elements = static_cast<size_t>(std::max<int64_t>(col.getNelem(), 1));

    appendColumnValues<int8_t>(col, values.re);
    appendColumnValues<uint8_t>(col, values.re);
    appendColumnValues<int16_t>(col, values.re);
    appendColumnValues<uint16_t>(col, values.re);
    appendColumnValues<int32_t>(col, values.re);
    appendColumnValues<uint32_t>(col, values.re);
    appendColumnValues<int64_t>(col, values.re);
    appendColumnValues<uint64_t>(col, values.re);
    appendColumnValues<float>(col, values.re);
    appendColumnValues<double>(col, values.re);
    appendColumnParts<float>(col, values);
    appendColumnParts<double>(col, values);

    if (values.re.size() != col.size() * values.elements)
        throw std::invalid_argument("Column '" + col.getName() + "' isn't numeric, or its rows differ in length.");
    return values;
}

template<typename... Ts>
static bool payloadIsAnyOf(const std::type_index& payload)
{
    return ((payload == typeid(Ts) || payload == typeid(std::vector<Ts>)) || ...);
}

// Integer, float, double or complex, scalar or vector.
static bool isNumericColumn(const FITSform& col)
{
    return payloadIsAnyOf<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t, int64_t, uint64_t,
                          float, double, FITSform::complex, FITSform::dblcomplex>(col.payloadType());
}

// What the latest plot of a table viewer used, to find the rows behind a click or a box.
struct PlotSource
{
    enum class Kind { None, Scatter, Density } kind = Kind::None;
    std::string xCol;             // Scatter: empty for the row index. Density: alone for one complex or 2-vector column.
    std::string yCol;
    bool        filtered = false;
    RowSet      rows;             // only if filtered
    size_t      nBinsX = 0;       // Density only, as the log flags
    size_t      nBinsY = 0;
    bool        logX = false;
    bool        logY = false;
};

// Calls fn(row, x, y) with every point of the plot, in the coordinates of its axes.
template<typename Fn>
static void forEachPlotPoint(const FITStable& tab, const PlotSource& source, Fn fn)
{
    ColumnValues x, y;
    if (!source.xCol.empty()) x = readColumnValues(*tab.getColumn(source.xCol));
    if (!source.yCol.empty()) y = readColumnValues(*tab.getColumn(source.yCol));

    // Same treatment as _QtDensityCore on log axes.
    const auto emitPoint = [&](size_t row, double px, double py) {
        if (source.logX) { if (px <= 0.0) return; px = std::log10(px); }
        if (source.logY) { if (py <= 0.0) return; py = std::log10(py); }
        fn(row, px, py);
    };

    const auto visitRow = [&](size_t row) {
        if (source.kind == PlotSource::Kind::Density) {
            if (source.yCol.empty()) {
                if (x.isComplex) emitPoint(row, x.re[row], x.im[row]);
                else             emitPoint(row, x.re[2 * row], x.re[2 * row + 1]);
            } else {
                emitPoint(row, x.re[row], y.re[row]);
            }
        } else if (source.xCol.empty()) {
            const size_t n = y.elements;
            for (size_t e = 0; e < n; ++e) {
                if (y.isComplex)  emitPoint(row, y.re[row * n + e], y.im[row * n + e]);
                else if (n == 1)  emitPoint(row, static_cast<double>(row + 1), y.re[row]);
                else              emitPoint(row, static_cast<double>(e), y.re[row * n + e]);
            }
        } else {
            const size_t n = y.elements;
            for (size_t e = 0; e < n; ++e)
                emitPoint(row, x.re[row * n + e], y.re[row * n + e]);
        }
    };

    if (source.filtered) {
        for (size_t row : source.rows.indices())
            visitRow(row);
    } else {
        const size_t nRows = tab.nrows();
        for (size_t row = 0; row < nRows; ++row)
            visitRow(row);
    }
}

using PlotPoints = std::vector<std::pair<double, double>>;

static RowSet rowsInBox(const FITStable& tab, const PlotSource& source, const QCPRange& xRange, const QCPRange& yRange, PlotPoints& marks)
{
    std::vector<size_t> rows;
    forEachPlotPoint(tab, source, [&](size_t row, double x, double y) {
        if (xRange.contains(x) && yRange.contains(y)) {
            rows.push_back(row);
            marks.emplace_back(x, y);
        }
    });
    return RowSet(std::move(rows));
}

// The point nearest to a position on screen, within a radius in pixels.
static RowSet rowNearest(const FITStable& tab, const PlotSource& source, QCPAxis* xAxis, QCPAxis* yAxis,
                         const QPoint& position, double radius, PlotPoints& marks)
{
    double best = radius * radius;
    bool   found = false;
    size_t bestRow = 0;
    std::pair<double, double> bestPoint;

    forEachPlotPoint(tab, source, [&](size_t row, double x, double y) {
        if (!std::isfinite(x) || !std::isfinite(y))
            return;
        const double dx = xAxis->coordToPixel(x) - position.x();
        const double dy = yAxis->coordToPixel(y) - position.y();
        const double distance = dx * dx + dy * dy;
        if (distance <= best) {
            best = distance;
            found = true;
            bestRow = row;
            bestPoint = {x, y};
        }
    });

    if (!found)
        return RowSet();
    marks.push_back(bestPoint);
    return RowSet(std::vector<size_t>{bestRow});
}

// The rows counted in the heatmap cell holding (cellX, cellY), binned as _QtDensityCore bins them.
static RowSet rowsInDensityCell(const FITStable& tab, const PlotSource& source, double cellX, double cellY)
{
    double xMin = std::numeric_limits<double>::infinity(), xMax = -xMin;
    double yMin = xMin, yMax = xMax;
    forEachPlotPoint(tab, source, [&](size_t, double x, double y) {
        if (std::isnan(x) || std::isnan(y))
            return;
        xMin = std::min(xMin, x);  xMax = std::max(xMax, x);
        yMin = std::min(yMin, y);  yMax = std::max(yMax, y);
    });
    if (xMin > xMax)
        return RowSet();

    const double xRange = (xMax > xMin) ? (xMax - xMin) : 1.0;
    const double yRange = (yMax > yMin) ? (yMax - yMin) : 1.0;

    // The very expression of _QtDensityCore, so that rounding puts a value in the same bin.
    const auto binOfX = [&](double x) { return std::floor(((x - xMin) / xRange) * static_cast<double>(source.nBinsX - 1) + 0.5); };
    const auto binOfY = [&](double y) { return std::floor(((y - yMin) / yRange) * static_cast<double>(source.nBinsY - 1) + 0.5); };

    const double wantedX = binOfX(cellX);
    const double wantedY = binOfY(cellY);
    if (wantedX < 0.0 || wantedX > static_cast<double>(source.nBinsX - 1) ||
        wantedY < 0.0 || wantedY > static_cast<double>(source.nBinsY - 1))
        return RowSet();

    std::vector<size_t> rows;
    forEachPlotPoint(tab, source, [&](size_t row, double x, double y) {
        if (binOfX(x) == wantedX && binOfY(y) == wantedY)
            rows.push_back(row);
    });
    return RowSet(std::move(rows));
}

// ---------------------------------------------------------------------------
// Catalogue overlay on a 2-D image
// ---------------------------------------------------------------------------

enum class CatalogueFrame { World, PixelFITS, PixelArray };

struct CatalogueChoice
{
    std::shared_ptr<const FITStable> table;
    std::string    xCol;
    std::string    yCol;   // empty: both coordinates come from xCol
    CatalogueFrame frame = CatalogueFrame::PixelFITS;
    std::string    filter;
    QString        label;  // names the catalogue to the user
};

struct CataloguePoint { size_t row; double x; double y; };

// Asks which table and columns to overlay, offering choice.filter as the row filter.
// The table is only read, into memory.
static bool askCatalogue(QWidget* parent, const std::string& fitsPath, bool imageHasWCS, CatalogueChoice& choice)
{
    QDialog dialog(parent);
    dialog.setWindowTitle("Catalogue overlay");
    auto* form = new QFormLayout(&dialog);

    auto* fileEdit = new QLineEdit(QString::fromStdString(fitsPath));
    fileEdit->setReadOnly(true);
    auto* browseBtn = new QPushButton("Browse...");
    auto* fileRow = new QHBoxLayout;
    fileRow->addWidget(fileEdit, 1);
    fileRow->addWidget(browseBtn);
    form->addRow("File:", fileRow);

    auto* hduCombo = new QComboBox;
    form->addRow("Table HDU:", hduCombo);

    auto* worldRadio = new QRadioButton("World coordinates, converted with the image's WCS");
    auto* fitsRadio  = new QRadioButton("Pixels, FITS (first pixel = 1)");
    auto* arrayRadio = new QRadioButton("Pixels, array (first pixel = 0)");
    worldRadio->setEnabled(imageHasWCS);
    if (!imageHasWCS)
        worldRadio->setToolTip("This image has no WCS calibration.");
    (imageHasWCS ? worldRadio : fitsRadio)->setChecked(true);
    auto* frameLayout = new QVBoxLayout;
    frameLayout->addWidget(worldRadio);
    frameLayout->addWidget(fitsRadio);
    frameLayout->addWidget(arrayRadio);
    form->addRow("Coordinates are:", frameLayout);

    auto* xCombo = new QComboBox;
    auto* yCombo = new QComboBox;
    yCombo->setToolTip("(from X column): a complex column gives real = X and imaginary = Y, "
                       "a 2-element vector column gives [0] = X and [1] = Y.");
    form->addRow("X column:", xCombo);
    form->addRow("Y column:", yCombo);

    auto* filterEdit = new QLineEdit(QString::fromStdString(choice.filter));
    filterEdit->setPlaceholderText("Optional row filter, e.g. phot_g_mean_mag < 18");
    form->addRow("Filter:", filterEdit);

    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    form->addRow(buttons);

    std::shared_ptr<const FITStable> table;
    std::vector<std::string> columnNames;

    const auto loadColumns = [&]() {
        table.reset();
        columnNames.clear();
        const QSignalBlocker blockX(xCombo), blockY(yCombo);
        xCombo->clear();
        yCombo->clear();
        yCombo->addItem("(from X column)");

        if (hduCombo->currentIndex() < 0)
            return;
        try {
            DataFile fm(fileEdit->text().toStdString());
            table = std::const_pointer_cast<const FITStable>(fm.GetTableAtIndex(hduCombo->currentData().toInt()));
        } catch (const std::exception& e) {
            QMessageBox::warning(&dialog, "Catalogue", QString::fromStdString(e.what()));
            return;
        }
        for (const auto& entry : table->listColumns()) {
            if (!isNumericColumn(*table->getColumn(entry[0])))
                continue;
            columnNames.push_back(entry[0]);
            const QString label = QString::fromStdString(entry[0])
                + (entry[2].empty() ? QString() : " [" + QString::fromStdString(entry[2]) + "]");
            xCombo->addItem(label);
            yCombo->addItem(label);
        }
        if (yCombo->count() > 2)
            yCombo->setCurrentIndex(2);
    };

    const auto loadHdus = [&]() {
        const QSignalBlocker block(hduCombo);
        hduCombo->clear();
        try {
            DataFile fm(fileEdit->text().toStdString());
            const std::vector<HduInfo> hdus = enumerateHDUs(fm);
            fm.Close();
            for (const auto& hdu : hdus)
                if (hdu.type == HduInfo::Type::BTBL || hdu.type == HduInfo::Type::TBL)
                    hduCombo->addItem(QString("HDU %1 - %2").arg(hdu.index).arg(QString::fromStdString(hdu.extname)), hdu.index);
        } catch (const std::exception& e) {
            QMessageBox::warning(&dialog, "Catalogue", QString::fromStdString(e.what()));
        }
        loadColumns();
    };

    QObject::connect(hduCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), &dialog, [&](int) { loadColumns(); });
    QObject::connect(browseBtn, &QPushButton::clicked, &dialog, [&]() {
        const QString path = QFileDialog::getOpenFileName(&dialog, "Catalogue file", fileEdit->text());
        if (!path.isEmpty()) {
            fileEdit->setText(path);
            loadHdus();
        }
    });
    QObject::connect(buttons, &QDialogButtonBox::accepted, &dialog, [&]() {
        if (!table || columnNames.empty())
            QMessageBox::warning(&dialog, "Catalogue", "This file has no table with a numeric column.");
        else
            dialog.accept();
    });
    QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);

    loadHdus();

    if (dialog.exec() != QDialog::Accepted)
        return false;

    choice.table  = table;
    choice.xCol   = columnNames[static_cast<size_t>(xCombo->currentIndex())];
    choice.yCol   = yCombo->currentIndex() > 0 ? columnNames[static_cast<size_t>(yCombo->currentIndex() - 1)] : std::string();
    choice.frame  = worldRadio->isChecked() ? CatalogueFrame::World
                  : fitsRadio->isChecked()  ? CatalogueFrame::PixelFITS
                                            : CatalogueFrame::PixelArray;
    choice.filter = filterEdit->text().trimmed().toStdString();
    choice.label  = QFileInfo(fileEdit->text()).fileName() + ", " + hduCombo->currentText();
    if (!choice.filter.empty())
        choice.label += " [" + filterEdit->text().trimmed() + "]";
    return true;
}

// The positions the catalogue's selected rows give, as read from its columns, NaN included.
static std::vector<CataloguePoint> readCataloguePositions(const CatalogueChoice& choice)
{
    const RowSet rows = choice.table->where(choice.filter);

    const FITSform& xCol = *choice.table->getColumn(choice.xCol);
    const ColumnValues x = readColumnValues(xCol);

    std::vector<CataloguePoint> points;

    if (choice.yCol.empty()) {
        if (!x.isComplex && x.elements != 2)
            throw std::invalid_argument("Column '" + choice.xCol + "' can't give both coordinates: "
                                        "only a complex or a 2-element vector column can. Choose a Y column.");
        for (size_t row : rows.indices()) {
            if (x.isComplex)
                for (size_t e = 0; e < x.elements; ++e)
                    points.push_back({row, x.re[row * x.elements + e], x.im[row * x.elements + e]});
            else
                points.push_back({row, x.re[2 * row], x.re[2 * row + 1]});
        }
        return points;
    }

    const ColumnValues y = readColumnValues(*choice.table->getColumn(choice.yCol));
    if (x.isComplex || y.isComplex)
        throw std::invalid_argument("A complex column holds both coordinates: choose it as X and leave Y on (from X column).");
    if (x.elements != y.elements)
        throw std::invalid_argument("Columns '" + choice.xCol + "' and '" + choice.yCol + "' don't have the same number of elements.");

    for (size_t row : rows.indices())
        for (size_t e = 0; e < x.elements; ++e)
            points.push_back({row, x.re[row * x.elements + e], y.re[row * x.elements + e]});
    return points;
}

// Turns positions into the 0-based pixels the image is drawn with; a position with no pixel is dropped.
static std::vector<CataloguePoint> cataloguePixels(const std::vector<CataloguePoint>& positions, CatalogueFrame frame, const FITScube& img)
{
    std::vector<CataloguePoint> pixels;
    pixels.reserve(positions.size());

    if (frame != CatalogueFrame::World) {
        const double shift = (frame == CatalogueFrame::PixelFITS) ? 1.0 : 0.0;
        for (const auto& p : positions)
            if (std::isfinite(p.x) && std::isfinite(p.y))
                pixels.push_back({p.row, p.x - shift, p.y - shift});
        return pixels;
    }

    std::vector<CataloguePoint> finite;
    worldVectors world;
    finite.reserve(positions.size());
    world.reserve(positions.size());
    for (const auto& p : positions) {
        if (std::isfinite(p.x) && std::isfinite(p.y)) {
            finite.push_back(p);
            world.push_back({p.x, p.y});
        }
    }
    if (world.empty())
        return pixels;

    // A GWCS gives NaN for a position whose pixel is outside its bounding box: dropped too
    const auto keep = [&](const CataloguePoint& from, const pixelCoords& p) {
        if (p.size() >= 2 && std::isfinite(p[0]) && std::isfinite(p[1]))
            pixels.push_back({from.row, p[0], p[1]});
    };

    try {
        const pixelVectors converted = img.World2PixelVector(world);
        for (size_t i = 0; i < finite.size(); ++i)
            keep(finite[i], converted[i]);
    } catch (const std::exception&) {
        // wcslib refuses the whole list for a single invalid position, such as one on the far
        // side of the sky for this projection: convert one by one and drop those.
        pixels.clear();
        for (size_t i = 0; i < finite.size(); ++i) {
            try {
                const pixelVectors converted = img.World2PixelVector(worldVectors{world[i]});
                keep(finite[i], converted[0]);
            } catch (const std::exception&) {
            }
        }
    }
    return pixels;
}

// World coordinates go to the WCS as they are, so say so if a column's unit isn't the axis's.
static bool confirmCatalogueUnits(QWidget* parent, const CatalogueChoice& choice, const FITScube& img)
{
    const auto units = img.WorldAxisUnits();
    const std::string yCol = choice.yCol.empty() ? choice.xCol : choice.yCol;
    const std::string xUnit = choice.table->getColumn(choice.xCol)->getUnit();
    const std::string yUnit = choice.table->getColumn(yCol)->getUnit();

    const auto differs = [&](const std::string& unit, size_t axis) {
        if (unit.empty() || axis >= units.size())
            return false;
        return QString::fromStdString(unit).trimmed() != QString::fromStdString(units[axis]).trimmed();
    };

    if (!differs(xUnit, 0) && !differs(yUnit, 1))
        return true;

    return QMessageBox::warning(parent, "Catalogue",
        QString("The columns' units (%1, %2) differ from the image's WCS axes (%3, %4). "
                "Values are not converted. Continue?")
            .arg(QString::fromStdString(xUnit), QString::fromStdString(yUnit),
                 QString::fromStdString(units.size() > 0 ? units[0] : std::string()),
                 QString::fromStdString(units.size() > 1 ? units[1] : std::string())),
        QMessageBox::Yes | QMessageBox::Cancel, QMessageBox::Cancel) == QMessageBox::Yes;
}

// Ten colours by ten markers: as many distinct styles as catalogues can be overlaid at once.
static constexpr size_t kMaxCatalogues = 100;

// The style of the index-th catalogue overlaid. Consecutive catalogues differ by both colour and
// marker, and each colour and marker pair comes once in kMaxCatalogues.
static QCPScatterStyle catalogueStyle(size_t index)
{
    static const QColor colours[] = {
        QColor(0, 255, 0),     QColor(255, 60, 60),   QColor(0, 230, 255),  QColor(255, 0, 255),
        QColor(255, 255, 0),   QColor(255, 150, 0),   QColor(90, 140, 255), QColor(255, 255, 255),
        QColor(255, 170, 200), QColor(170, 110, 255)
    };
    static const QCPScatterStyle::ScatterShape markers[] = {
        QCPScatterStyle::ssCircle,    QCPScatterStyle::ssSquare,
        QCPScatterStyle::ssDiamond,   QCPScatterStyle::ssTriangle,
        QCPScatterStyle::ssTriangleInverted, QCPScatterStyle::ssCross,
        QCPScatterStyle::ssPlus,      QCPScatterStyle::ssStar,
        QCPScatterStyle::ssCrossCircle, QCPScatterStyle::ssPlusSquare
    };
    const size_t colour = index % 10;
    const size_t marker = (colour + index / 10) % 10;
    return QCPScatterStyle(markers[marker], QPen(colours[colour], 1.2), Qt::NoBrush, 9);
}

// The marker alone, on a dark background so that the light colours show in a menu.
static QIcon catalogueIcon(const QCPScatterStyle& style)
{
    QPixmap pixmap(16, 16);
    pixmap.fill(QColor(40, 40, 40));
    QCPPainter painter(&pixmap);
    painter.setAntialiasing(true);
    style.applyTo(&painter, QPen());
    style.drawShape(&painter, 8.0, 8.0);
    painter.end();
    return QIcon(pixmap);
}

// Adds to the image's toolbar what overlays catalogues, and clicking a marker to see its row.
static void addCatalogueTools(QMainWindow* win, QCustomPlot* plot,
                              const std::shared_ptr<const FITScube>& img, const std::string& fitsPath)
{
    auto* toolbar = win->findChild<QToolBar*>();
    if (!toolbar)
        return;

    auto* catalogueBtn = new QPushButton("Catalogue...", toolbar);
    catalogueBtn->setFlat(true);
    auto* keepCheck = new QCheckBox("Keep previous", toolbar);
    keepCheck->setToolTip(QString("Overlay the next catalogue over those displayed, with its own "
                                  "colour and marker (up to %1 catalogues).").arg(kMaxCatalogues));
    auto* clearBtn = new QPushButton("Clear catalogues", toolbar);
    clearBtn->setFlat(true);
    clearBtn->setEnabled(false);
    auto* keyBtn = new QToolButton(toolbar);
    keyBtn->setText("Catalogues");
    keyBtn->setToolTip("Which colour and marker each catalogue has");
    keyBtn->setPopupMode(QToolButton::InstantPopup);
    keyBtn->setEnabled(false);
    auto* keyMenu = new QMenu(keyBtn);
    keyBtn->setMenu(keyMenu);
    auto* countLabel = new QLabel(toolbar);

    toolbar->addSeparator();
    toolbar->addWidget(catalogueBtn);
    toolbar->addWidget(keepCheck);
    toolbar->addWidget(clearBtn);
    toolbar->addWidget(keyBtn);
    toolbar->addWidget(countLabel);

    struct Overlay
    {
        QCPGraph* graph = nullptr;
        std::shared_ptr<const FITStable> table;
        std::vector<CataloguePoint> points;
        QString label;
    };
    auto overlays = std::make_shared<std::vector<Overlay>>();
    auto lastFilter = std::make_shared<std::string>();
    QCPAxisRect* image = plot->axisRect(0);

    const auto clearOverlays = [=]() {
        for (const auto& overlay : *overlays)
            plot->removeGraph(overlay.graph);
        overlays->clear();
        keyMenu->clear();
        keyBtn->setText("Catalogues");
        keyBtn->setEnabled(false);
        countLabel->clear();
        clearBtn->setEnabled(false);
        plot->replot();
    };
    QObject::connect(clearBtn, &QPushButton::clicked, win, [=]() { clearOverlays(); });

    QObject::connect(catalogueBtn, &QPushButton::clicked, win, [=]() {
        if (keepCheck->isChecked() && overlays->size() >= kMaxCatalogues) {
            QMessageBox::warning(win, "Catalogue",
                QString("%1 catalogues are already displayed: clear them, or untick Keep previous.").arg(kMaxCatalogues));
            return;
        }

        CatalogueChoice choice;
        choice.filter = *lastFilter;
        if (!askCatalogue(win, fitsPath, img->HasWorldCoordinates(), choice))
            return;
        *lastFilter = choice.filter;

        std::vector<CataloguePoint> positions, pixels;
        try {
            if (choice.frame == CatalogueFrame::World && !confirmCatalogueUnits(win, choice, *img))
                return;
            positions = readCataloguePositions(choice);
            pixels    = cataloguePixels(positions, choice.frame, *img);
        } catch (const std::exception& e) {
            QMessageBox::warning(win, "Catalogue", QString::fromStdString(e.what()));
            return;
        }

        const double nx = static_cast<double>(img->Size(1));
        const double ny = static_cast<double>(img->Size(2));
        pixels.erase(std::remove_if(pixels.begin(), pixels.end(), [=](const CataloguePoint& p) {
            return !(p.x >= -0.5 && p.x < nx - 0.5 && p.y >= -0.5 && p.y < ny - 0.5);
        }), pixels.end());

        if (!keepCheck->isChecked())
            clearOverlays();

        QVector<double> xs, ys;
        xs.reserve(static_cast<int>(pixels.size()));
        ys.reserve(static_cast<int>(pixels.size()));
        for (const auto& p : pixels) {
            xs.push_back(p.x);
            ys.push_back(p.y);
        }

        const QCPScatterStyle style = catalogueStyle(overlays->size());

        Overlay overlay;
        overlay.graph = plot->addGraph(image->axis(QCPAxis::atBottom), image->axis(QCPAxis::atLeft));
        overlay.graph->setLineStyle(QCPGraph::lsNone);
        overlay.graph->setScatterStyle(style);
        overlay.graph->setData(xs, ys);
        overlay.table  = choice.table;
        overlay.points = std::move(pixels);
        overlay.label  = choice.label;

        const QString sources = QString("%1 of %2 sources in image").arg(overlay.points.size()).arg(positions.size());
        QString keyText = overlay.label + " - " + sources;
        keyMenu->addAction(catalogueIcon(style), keyText.replace("&", "&&"));
        overlays->push_back(std::move(overlay));

        keyBtn->setText(QString("Catalogues (%1)").arg(overlays->size()));
        keyBtn->setEnabled(true);
        countLabel->setText(" " + sources);
        clearBtn->setEnabled(true);
        plot->replot();
    });

    // A click is a press and release at the same place; a drag pans the image.
    auto pressedAt = std::make_shared<QPoint>();
    QObject::connect(plot, &QCustomPlot::mousePress, win, [=](QMouseEvent*) {
        *pressedAt = plot->mapFromGlobal(QCursor::pos());
    });
    QObject::connect(plot, &QCustomPlot::mouseRelease, win, [=](QMouseEvent* event) {
        const QPoint position = plot->mapFromGlobal(QCursor::pos());
        if (event->button() != Qt::LeftButton || overlays->empty() ||
            (position - *pressedAt).manhattanLength() > 3 || plot->axisRectAt(position) != image)
            return;

        QCPAxis* xAxis = image->axis(QCPAxis::atBottom);
        QCPAxis* yAxis = image->axis(QCPAxis::atLeft);

        double best = 10.0 * 10.0;
        const Overlay* nearestOverlay = nullptr;
        const CataloguePoint* nearest = nullptr;
        for (const auto& overlay : *overlays) {
            for (const auto& p : overlay.points) {
                const double dx = xAxis->coordToPixel(p.x) - position.x();
                const double dy = yAxis->coordToPixel(p.y) - position.y();
                const double distance = dx * dx + dy * dy;
                if (distance <= best) {
                    best = distance;
                    nearestOverlay = &overlay;
                    nearest = &p;
                }
            }
        }
        if (nearest)
            showSelectedRows(win, nearestOverlay->table, RowSet(std::vector<size_t>{nearest->row}),
                             nearestOverlay->label + QString(", row %1").arg(nearest->row + 1));
    });
}

static QWidget* makeFITSTableViewer(const HduInfo& hdu, const std::string& fitsPath)
{
    // --- Load the table ---
    DataFile fm(fitsPath);
    // A copy in memory: the file is opened read-only and nothing here writes to it. Only
    // computeColumn() goes through the non-const pointer, and it stays in this copy too.
    auto table = fm.GetTableAtIndex(hdu.index);
    auto tab = std::const_pointer_cast<const FITStable>(table);

    // --- Column metadata ---
    const FITStable::clist cols = tab->listColumns(); // each entry: {name, TFORM, unit}
    const size_t nCols = tab->ncols();
    const size_t nRows = tab->nrows();

    // Build separate plottable lists for X and Y combo boxes.
    // plottableX: real/integer scalar columns (no complex) — used for Scatter X, Density Z weight, Histogram.
    // plottableY: all numeric columns including complex — used for Scatter Y, Density X/Y.
    // Complex columns are allowed in Density X/Y selectors; when one axis is left at "-- (none) --"
    // and the other is complex, QtDensityPlot's single-column overload splits real/imag parts.
    // Shared with the plot lambdas, which must see the columns computed later on.
    struct PlotCol { std::string name; std::string unit; };
    auto plottableX = std::make_shared<std::vector<PlotCol>>();
    auto plottableY = std::make_shared<std::vector<PlotCol>>();
    auto categoryCols = std::make_shared<std::vector<PlotCol>>();
    auto vectorYCols = std::make_shared<std::vector<PlotCol>>();
    auto vectorXCols = std::make_shared<std::vector<PlotCol>>();

    for (const auto& entry : cols) {
        const std::string& name = entry[0];
        const std::string& unit = entry[2];
        const auto& col = tab->getColumn(name);
        if (isPlottableAsX(col)) plottableX->push_back({name, unit});
        if (isPlottableAsY(col)) plottableY->push_back({name, unit});
        if (isCategoryColumn(*col)) categoryCols->push_back({name, unit});
        if (col->getNelem() > 1 && isPlottableAsY(col)) vectorYCols->push_back({name, unit});
        if (col->getNelem() > 1 && isPlottableAsX(col)) vectorXCols->push_back({name, unit});
    }

    // Label helper: "name [unit]" or just "name" when no unit
    auto colLabel = [](const PlotCol& c) -> QString {
        return c.unit.empty()
            ? QString::fromStdString(c.name)
            : QString::fromStdString(c.name) + " [" + QString::fromStdString(c.unit) + "]";
    };

    // --- Window shell ---
    auto* win = new QMainWindow;
    win->setWindowTitle(QString("HDU %1 - %2 (%3)")
        .arg(hdu.index)
        .arg(QString::fromStdString(hdu.extname))
        .arg(hdu.type == HduInfo::Type::BTBL ? "Binary Table" : "ASCII Table"));

    // --- Left panel: table summary (unchanged) ---
    auto* leftPanel  = new QWidget;
    auto* leftLayout = new QVBoxLayout(leftPanel);
    auto* columnCountLabel = new QLabel(QString("Columns: %1").arg(nCols));
    leftLayout->addWidget(columnCountLabel);
    leftLayout->addWidget(new QLabel(QString("Rows: %1").arg(nRows)));

    auto* colTable = new QTableWidget(static_cast<int>(nCols), 2);
    colTable->setHorizontalHeaderLabels({"Name", "Type"});
    colTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    colTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    colTable->verticalHeader()->setVisible(false);
    for (int i = 0; i < static_cast<int>(nCols); ++i) {
        const auto& c = cols[static_cast<size_t>(i)]; // c[0]=name, c[1]=TFORM, c[2]=unit
        colTable->setItem(i, 0, new QTableWidgetItem(QString::fromStdString(c[0])));
        QString typeStr = QString::fromStdString(c[1]);
        if (!c[2].empty())
            typeStr += " [" + QString::fromStdString(c[2]) + "]";
        colTable->setItem(i, 1, new QTableWidgetItem(typeStr));
    }
    enableCellCopy(colTable);
    leftLayout->addWidget(colTable, 1);
    leftPanel->setMinimumWidth(120);

    // --- Centre panel: plot area ---
    auto* customPlot = new QCustomPlot;
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes);
    customPlot->axisRect()->setupFullAxesBox(true);
    customPlot->legend->setVisible(false);  // axis labels identify columns; legend causes layout crash on first show
    customPlot->setMinimumWidth(200);
    // A highlighted category is drawn over the other plots, and the selection rings over it.
    customPlot->addLayer("categoryHighlight", customPlot->layer("main"), QCustomPlot::limAbove);
    customPlot->addLayer("selectionMarks", customPlot->layer("categoryHighlight"), QCustomPlot::limAbove);

    // Row filter above the plot, with the count of rows displayed beside it
    auto* filterEdit = new QLineEdit;
    filterEdit->setPlaceholderText("Row filter, e.g. (OBJNAME=~WD* || SCSID>0) && RA >= 10.5");
    filterEdit->setClearButtonEnabled(true);
    auto* rowCountLabel = new QLabel(QString("0 / %1 rows").arg(nRows));
    rowCountLabel->setToolTip("Rows displayed by the latest plot / rows in the table");

    auto* filterLayout = new QHBoxLayout;
    filterLayout->addWidget(new QLabel("Filter:"));
    filterLayout->addWidget(filterEdit, 1);
    filterLayout->addWidget(rowCountLabel);
    auto* selectionLabel = new QLabel;
    filterLayout->addWidget(selectionLabel);

    // New column computed from the others: kept in memory until the table is reopened
    auto* newColumnEdit = new QLineEdit;
    newColumnEdit->setPlaceholderText("In memory only, e.g. COLOR = MagG - MagPb");
    newColumnEdit->setClearButtonEnabled(true);
    newColumnEdit->setToolTip("NAME = expression, computed on every row and kept until the table is reopened; "
                              "the file is never changed.\n"
                              "Operators: + - * / ^ ( )\n"
                              "Functions: sqrt exp log log10 sin cos tan asin acos atan abs pow(x, y)\n"
                              "Complex columns: real img abs arg\n"
                              "Vector columns work element by element.");
    auto* newColumnBtn = new QPushButton("Add");

    auto* newColumnLayout = new QHBoxLayout;
    newColumnLayout->addWidget(new QLabel("New column:"));
    newColumnLayout->addWidget(newColumnEdit, 1);
    newColumnLayout->addWidget(newColumnBtn);

    auto* centerPanel  = new QWidget;
    auto* centerLayout = new QVBoxLayout(centerPanel);
    centerLayout->setContentsMargins(0, 0, 0, 0);
    centerLayout->addLayout(filterLayout);
    centerLayout->addLayout(newColumnLayout);
    centerLayout->addWidget(customPlot, 1);

    // --- Right panel: QTabWidget with three tabs ---
    auto* tabWidget = new QTabWidget;
    tabWidget->setMinimumWidth(160);

    // ---- Tab 0: Scatter ----
    auto* scatterTab    = new QWidget;
    auto* scatterLayout = new QVBoxLayout(scatterTab);
    auto* scatterForm   = new QFormLayout;
    auto* xScatterCombo = new QComboBox;
    auto* yScatterCombo = new QComboBox;

    // X: "(Row index)" sentinel first, then real/integer columns (no complex)
    xScatterCombo->addItem("(Row index)");
    for (const auto& c : *plottableX)
        xScatterCombo->addItem(colLabel(c));

    // Y: all numeric columns including complex; no sentinel — column required
    for (const auto& c : *plottableY)
        yScatterCombo->addItem(colLabel(c));

    scatterForm->addRow("X axis:", xScatterCombo);
    scatterForm->addRow("Y axis:", yScatterCombo);
    scatterLayout->addLayout(scatterForm);
    scatterLayout->addStretch();
    tabWidget->addTab(scatterTab, "Scatter");

    // ---- Tab 1: Density Plot ----
    auto* densTab    = new QWidget;
    auto* densLayout = new QVBoxLayout(densTab);
    auto* densForm   = new QFormLayout;
    auto* xDensCombo = new QComboBox;
    auto* yDensCombo = new QComboBox;
    auto* zDensCombo = new QComboBox;

    // X/Y: "-- (none) --" sentinel first (index 0), then all numeric columns incl. complex
    xDensCombo->addItem("-- (none) --");
    yDensCombo->addItem("-- (none) --");
    for (const auto& c : *plottableY) {
        xDensCombo->addItem(colLabel(c));
        yDensCombo->addItem(colLabel(c));
    }

    // Z/weight: "-- (none) --" sentinel first, then real/integer columns only
    zDensCombo->addItem("-- (none) --");
    for (const auto& c : *plottableX)
        zDensCombo->addItem(colLabel(c));

    auto* densNBinsX = new QSpinBox;
    auto* densNBinsY = new QSpinBox;
    densNBinsX->setRange(2, 2000);  densNBinsX->setValue(100);
    densNBinsY->setRange(2, 2000);  densNBinsY->setValue(100);

    auto* densLogX = new QCheckBox("Log X");
    auto* densLogY = new QCheckBox("Log Y");
    auto* densLogZ = new QCheckBox("Log color scale");
    auto* densLogRow = new QHBoxLayout;
    densLogRow->addWidget(densLogX);
    densLogRow->addWidget(densLogY);
    densLogRow->addWidget(densLogZ);

    densForm->addRow("X axis:", xDensCombo);
    densForm->addRow("Y axis:", yDensCombo);
    densForm->addRow("Weight (optional):", zDensCombo);
    densForm->addRow("Bins X:", densNBinsX);
    densForm->addRow("Bins Y:", densNBinsY);
    densLayout->addLayout(densForm);
    densLayout->addLayout(densLogRow);
    densLayout->addStretch();
    tabWidget->addTab(densTab, "Density");

    // ---- Tab 2: Histogram ----
    auto* histTab    = new QWidget;
    auto* histLayout = new QVBoxLayout(histTab);
    auto* histForm   = new QFormLayout;
    auto* histCombo  = new QComboBox;

    // Column: real/integer only; no sentinel — column required
    for (const auto& c : *plottableX)
        histCombo->addItem(colLabel(c));

    auto* histNBins = new QSpinBox;
    histNBins->setRange(0, 10000);
    histNBins->setValue(0);           // 0 → auto (Sturges rule)
    histNBins->setSpecialValueText("auto");

    auto* histLogX  = new QCheckBox("Log bins (X)");
    auto* histLogY  = new QCheckBox("Log Y axis");
    auto* histLogRow = new QHBoxLayout;
    histLogRow->addWidget(histLogX);
    histLogRow->addWidget(histLogY);

    histForm->addRow("Column:", histCombo);
    histForm->addRow("Bins:", histNBins);
    histLayout->addLayout(histForm);
    histLayout->addLayout(histLogRow);
    histLayout->addStretch();
    tabWidget->addTab(histTab, "Histogram");

    // ---- Tab 3: Vector row ----
    auto* vectorTab    = new QWidget;
    auto* vectorLayout = new QVBoxLayout(vectorTab);
    auto* vectorForm   = new QFormLayout;
    auto* vecYCombo    = new QComboBox;
    auto* vecXCombo    = new QComboBox;

    for (const auto& c : *vectorYCols)
        vecYCombo->addItem(colLabel(c));

    vecXCombo->addItem("(Element index)");
    for (const auto& c : *vectorXCols)
        vecXCombo->addItem(colLabel(c));

    auto* vecRowSpin = new QSpinBox;
    vecRowSpin->setRange(1, static_cast<int>(std::clamp<size_t>(nRows, 1, std::numeric_limits<int>::max())));
    vecRowSpin->setKeyboardTracking(false); // a row typed digit by digit is only drawn once
    auto* vecPrevBtn = new QPushButton("Previous");
    auto* vecNextBtn = new QPushButton("Next");
    auto* vecStepRow = new QHBoxLayout;
    vecStepRow->addWidget(vecPrevBtn);
    vecStepRow->addWidget(vecNextBtn);

    auto* vecLineBox = new QCheckBox("Overlay line");

    vectorForm->addRow("X:", vecXCombo);
    vectorForm->addRow("Y:", vecYCombo);
    vectorForm->addRow("Row:", vecRowSpin);
    vectorLayout->addLayout(vectorForm);
    vectorLayout->addLayout(vecStepRow);
    vectorLayout->addWidget(vecLineBox);
    vectorLayout->addStretch();
    tabWidget->addTab(vectorTab, "Vector row");
    tabWidget->setTabEnabled(3, !vectorYCols->empty());

    // ---- Category: splits Scatter and Histogram into one colored series per category ----
    auto* categoryCombo = new QComboBox;
    categoryCombo->addItem("-- (none) --");
    for (const auto& c : *categoryCols)
        categoryCombo->addItem(colLabel(c));

    auto* epsilonEdit = new QLineEdit("0");
    auto* epsilonValidator = new QDoubleValidator(0.0, std::numeric_limits<double>::max(), 1000, epsilonEdit);
    epsilonValidator->setNotation(QDoubleValidator::ScientificNotation);
    epsilonValidator->setLocale(QLocale::c());
    epsilonEdit->setValidator(epsilonValidator);
    epsilonEdit->setEnabled(false);
    epsilonEdit->setToolTip("Float and double columns: each category spans value ± epsilon, "
                            "with the next one starting right above. 0 gives one category per distinct value.");

    // Plain list rather than customPlot->legend, which crashes on first show (see above).
    auto* categoryKey = new QListWidget;
    categoryKey->setSelectionMode(QAbstractItemView::NoSelection);
    categoryKey->setToolTip("Click a category to highlight it on the plot, and again to show them all.");
    categoryKey->setVisible(false);

    auto* categoryBox    = new QGroupBox("Category");
    auto* categoryLayout = new QVBoxLayout(categoryBox);
    auto* categoryForm   = new QFormLayout;
    categoryForm->addRow("Column:", categoryCombo);
    categoryForm->addRow("Epsilon:", epsilonEdit);
    categoryLayout->addLayout(categoryForm);
    categoryLayout->addWidget(categoryKey, 1);

    // --- Overlay / layering state ---
    // hasContent tracks whether the canvas currently has any plotted data.
    // contourColorIdx and scatterColorIdx cycle through distinct colors on each
    // successive Plot call in Retain mode so overlaid plots are visually distinct.
    auto hasContent      = std::make_shared<bool>(false);
    auto contourColorIdx = std::make_shared<int>(0);
    auto scatterColorIdx = std::make_shared<int>(0);

    // Color palettes cycled in Retain mode
    static const QColor scatterPalette[] = {
        Qt::blue, Qt::red, Qt::darkGreen, Qt::magenta, Qt::cyan, QColor(255,128,0)
    };
    static const QColor contourPalette[] = {
        Qt::white, Qt::yellow, Qt::cyan, Qt::magenta, Qt::red, Qt::green
    };
    static const int nScatterColors = 6;
    static const int nContourColors = 6;

    auto showRowCount = [=](size_t displayed) {
        rowCountLabel->setText(QString("%1 / %2 rows").arg(displayed).arg(nRows));
    };

    // Called once a plot has made it to the canvas; rows is null when no filter applies.
    auto markPlotted = [=](const RowSet* rows) {
        *hasContent = true;
        showRowCount(rows ? rows->size() : nRows);
    };

    // The graphs drawn for each category of the latest plot, in the order of its list.
    // QPointer: the plot owns them and may delete them.
    using CategoryGraphs = std::vector<QPointer<QCPGraph>>;
    struct KeyedCategory { QColor color; CategoryGraphs graphs; };
    auto keyedCategories     = std::make_shared<std::vector<KeyedCategory>>();
    auto highlightedCategory = std::make_shared<int>(-1);

    // Lists the categories of the latest plot with their color; hidden when it has none.
    auto showCategoryKey = [=](const std::vector<TableCategory>& categories, const std::vector<CategoryGraphs>& graphs) {
        categoryKey->clear();
        keyedCategories->clear();
        *highlightedCategory = -1;
        for (size_t i = 0; i < categories.size(); ++i) {
            QPixmap swatch(12, 12);
            swatch.fill(categories[i].color);
            categoryKey->addItem(new QListWidgetItem(QIcon(swatch),
                QString("%1 (%2 rows)").arg(categories[i].label, QString::number(categories[i].rows.size()))));
            keyedCategories->push_back({categories[i].color, graphs[i]});
        }
        categoryKey->setVisible(!categories.empty());
    };

    // Draws one category in its color above the others, dimmed; -1 gives all of them their color back.
    auto highlightCategory = [=](int category) {
        if (category == *highlightedCategory)
            return;
        *highlightedCategory = category;

        const QColor dimmed(215, 215, 215);
        for (int i = 0; i < static_cast<int>(keyedCategories->size()); ++i) {
            const KeyedCategory& keyed = (*keyedCategories)[static_cast<size_t>(i)];
            for (const auto& graph : keyed.graphs) {
                if (!graph)
                    continue;
                graph->setPen(QPen(category < 0 || i == category ? keyed.color : dimmed));
                // In list order, which restores the drawing order after a highlight.
                graph->setLayer(i == category ? "categoryHighlight" : "main");
            }
            QListWidgetItem* item = categoryKey->item(i);
            QFont font = item->font();
            font.setBold(i == category);
            item->setFont(font);
        }
        customPlot->replot(QCustomPlot::rpQueuedReplot);
    };

    QObject::connect(categoryKey, &QListWidget::itemClicked, win, [=](QListWidgetItem* item) {
        const int category = categoryKey->row(item);
        highlightCategory(category == *highlightedCategory ? -1 : category);
    });

    // What the latest plot used, to find the rows behind a click or a box, and the
    // ring marking the points selected on it. No plotted value is kept here.
    auto plotted        = std::make_shared<PlotSource>();
    auto selectionMarks = std::make_shared<QCPGraph*>(nullptr);

    auto clearSelection = [=]() {
        if (*selectionMarks)
            customPlot->removeGraph(*selectionMarks);
        *selectionMarks = nullptr;
        selectionLabel->clear();
    };

    // Remembers a plot just drawn from these columns, for the selection.
    auto rememberPlot = [=](PlotSource::Kind kind, const std::string& xCol, const std::string& yCol, const RowSet* rows) {
        PlotSource source;
        source.kind = kind;
        source.xCol = xCol;
        source.yCol = yCol;
        if (rows) {
            source.filtered = true;
            source.rows = *rows;
        }
        *plotted = std::move(source);
    };

    // --- clearCanvas lambda ---
    // Removes the QCPColorScale (column 1 of plotLayout), all plottables, and all
    // annotation items (contour labels).  Resets the color cycling counters.
    auto clearCanvas = [=]() {
        if (customPlot->plotLayout()->hasElement(0, 1)) {
            customPlot->plotLayout()->remove(customPlot->plotLayout()->element(0, 1));
            customPlot->plotLayout()->simplify();
        }
        *selectionMarks = nullptr; // about to be deleted with the other plottables
        selectionLabel->clear();
        *plotted = PlotSource();
        customPlot->clearPlottables();
        customPlot->clearItems();
        *hasContent      = false;
        *contourColorIdx = 0;
        *scatterColorIdx = 0;
        showRowCount(0);
        showCategoryKey({}, {});
        customPlot->replot();
    };

    // --- doScatter lambda ---
    auto doScatter = [=](const RowSet* rows, const std::vector<TableCategory>& categories) {
        const int xIdx = xScatterCombo->currentIndex();
        const bool useRowIdx = (xIdx == 0);
        const std::string xKey = useRowIdx
            ? std::string()
            : (*plottableX)[static_cast<size_t>(xIdx - 1)].name;
        const std::string yKey =
            (*plottableY)[static_cast<size_t>(yScatterCombo->currentIndex())].name;

        std::vector<std::shared_ptr<QCPGraph>> graphs;
        std::vector<CategoryGraphs> categoryGraphs(categories.size());
        try {
            if (!categories.empty()) {
                for (size_t i = 0; i < categories.size(); ++i) {
                    const RowSet& catRows = categories[i].rows;
                    auto catGraphs = useRowIdx
                        ? QtFITSviewer::QtScatter(customPlot, tab, yKey, catRows)
                        : QtFITSviewer::QtScatter(customPlot, tab, xKey, yKey, catRows);
                    for (auto& g : catGraphs) {
                        g->setPen(QPen(categories[i].color));
                        categoryGraphs[i].push_back(g.get());
                        graphs.push_back(g);
                    }
                }
            } else if (rows)
                graphs = useRowIdx
                    ? QtFITSviewer::QtScatter(customPlot, tab, yKey, *rows)
                    : QtFITSviewer::QtScatter(customPlot, tab, xKey, yKey, *rows);
            else
                graphs = useRowIdx
                    ? QtFITSviewer::QtScatter(customPlot, tab, yKey)
                    : QtFITSviewer::QtScatter(customPlot, tab, xKey, yKey);
        } catch (const std::exception& e) {
            QMessageBox::critical(nullptr, "Scatter error", QString::fromStdString(e.what()));
            return;
        }

        const QColor col = scatterPalette[(*scatterColorIdx)++ % nScatterColors];
        for (auto& g : graphs) {
            g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 4));
            if (categories.empty())
                g->setPen(QPen(col));
        }

        const QString xLabel = useRowIdx
            ? QString("Row index")
            : colLabel((*plottableX)[static_cast<size_t>(xIdx - 1)]);
        const QString yLabel = colLabel((*plottableY)[static_cast<size_t>(yScatterCombo->currentIndex())]);
        customPlot->xAxis->setLabel(xLabel);
        customPlot->yAxis->setLabel(yLabel);

        customPlot->rescaleAxes();
        customPlot->replot(QCustomPlot::rpQueuedReplot);
        markPlotted(rows);
        showCategoryKey(categories, categoryGraphs);
        rememberPlot(PlotSource::Kind::Scatter, xKey, yKey, rows);
    };

    // --- doHistogram lambda ---
    auto doHistogram = [=](const RowSet* rows, const std::vector<TableCategory>& categories) {
        if (plottableX->empty()) {
            QMessageBox::warning(nullptr, "Histogram", "No suitable columns available.");
            return;
        }
        const std::string colName =
            (*plottableX)[static_cast<size_t>(histCombo->currentIndex())].name;

        const size_t nBins   = static_cast<size_t>(histNBins->value()); // 0 = auto
        const bool   logBins = histLogX->isChecked();
        const bool   logY    = histLogY->isChecked();

        std::vector<std::shared_ptr<QCPGraph>> graphs;
        try {
            if (!categories.empty()) {
                std::vector<RowSet> groups;
                groups.reserve(categories.size());
                for (const auto& category : categories)
                    groups.push_back(category.rows);
                graphs = QtFITSviewer::QtHistogram(customPlot, tab, colName, groups, nBins);
            } else {
                graphs = rows
                    ? QtFITSviewer::QtHistogram(customPlot, tab, colName, *rows, nBins)
                    : QtFITSviewer::QtHistogram(customPlot, tab, colName, nBins);
            }
        } catch (const std::exception& e) {
            QMessageBox::critical(nullptr, "Histogram error", QString::fromStdString(e.what()));
            return;
        }

        // Apply log scale types after the fact (QtHistogram's FITStable overload
        // hardcodes logFlag=0; we set the axis scale directly on the plot).
        customPlot->xAxis->setScaleType(logBins ? QCPAxis::stLogarithmic : QCPAxis::stLinear);
        customPlot->yAxis->setScaleType(logY    ? QCPAxis::stLogarithmic : QCPAxis::stLinear);

        const QColor col = scatterPalette[(*scatterColorIdx)++ % nScatterColors];
        std::vector<CategoryGraphs> categoryGraphs(categories.size());
        for (size_t i = 0; i < graphs.size(); ++i) {
            graphs[i]->setPen(QPen(categories.empty() ? col : categories[i].color));
            if (!categories.empty())
                categoryGraphs[i].push_back(graphs[i].get());
        }

        customPlot->xAxis->setLabel(colLabel((*plottableX)[static_cast<size_t>(histCombo->currentIndex())]));
        customPlot->yAxis->setLabel("count");

        customPlot->rescaleAxes();
        customPlot->replot(QCustomPlot::rpQueuedReplot);
        markPlotted(rows);
        showCategoryKey(categories, categoryGraphs);
    };

    // --- doDensity lambda ---
    // When Retain is unchecked: draws a filled QCPColorMap (QtDensityPlot).
    // When Retain is checked: computes the density grid, then discards the color map
    // and overlays QCustomContours contour lines instead, using a new uniform color
    // on each call so successive overlaid contour layers are visually distinct.
    auto doDensity = [=](bool retain, const RowSet* rows) {
        const int xDensIdx = xDensCombo->currentIndex(); // 0 = "-- (none) --"
        const int yDensIdx = yDensCombo->currentIndex();
        const bool xIsNone = (xDensIdx == 0);
        const bool yIsNone = (yDensIdx == 0);

        // Validate: at least one axis must be selected
        if (xIsNone && yIsNone) {
            QMessageBox::warning(nullptr, "Density Plot",
                "Select at least one axis column.");
            return;
        }
        // If one axis is "-- (none) --" the other must be complex so the
        // single-column overload of QtDensityPlot can split real/imag parts.
        if (xIsNone || yIsNone) {
            const std::string checkCol = xIsNone
                ? (*plottableY)[static_cast<size_t>(yDensIdx - 1)].name
                : (*plottableY)[static_cast<size_t>(xDensIdx - 1)].name;
            const DSL::dtype t = tab->getColumn(checkCol)->getType();
            if (t != DSL::dtype::tcplx && t != DSL::dtype::tdbcplx) {
                QMessageBox::warning(nullptr, "Density Plot",
                    "A single-axis selection requires a complex-type column "
                    "(real and imaginary parts are used as X and Y).");
                return;
            }
        }

        // Compute the density map (always needed: either to display directly, or
        // to extract the binned grid for QCustomContours).
        const size_t  nBinsX  = static_cast<size_t>(densNBinsX->value());
        const size_t  nBinsY  = static_cast<size_t>(densNBinsY->value());
        // logFlag bitmask: 0x1=log₁₀(x), 0x2=log₁₀(y), 0x4=log₁₀(bin count/z)
        const uint8_t logFlag = static_cast<uint8_t>(
            (densLogX->isChecked() ? 0x1 : 0) |
            (densLogY->isChecked() ? 0x2 : 0) |
            (densLogZ->isChecked() ? 0x4 : 0));

        std::shared_ptr<QCPColorMap> colorMap;
        try {
            if (xIsNone || yIsNone) {
                // Single complex column: real→X, imag→Y
                const std::string singleCol = xIsNone
                    ? (*plottableY)[static_cast<size_t>(yDensIdx - 1)].name
                    : (*plottableY)[static_cast<size_t>(xDensIdx - 1)].name;
                colorMap = rows
                    ? QtFITSviewer::QtDensityPlot(customPlot, tab, singleCol, *rows,
                                                  nBinsX, nBinsY, logFlag)
                    : QtFITSviewer::QtDensityPlot(customPlot, tab, singleCol,
                                                  nBinsX, nBinsY, logFlag);
            } else {
                const std::string xCol = (*plottableY)[static_cast<size_t>(xDensIdx - 1)].name;
                const std::string yCol = (*plottableY)[static_cast<size_t>(yDensIdx - 1)].name;
                const int zIdx = zDensCombo->currentIndex();
                const std::string zCol = (zIdx == 0)
                    ? std::string()
                    : (*plottableX)[static_cast<size_t>(zIdx - 1)].name;
                colorMap = rows
                    ? QtFITSviewer::QtDensityPlot(customPlot, tab, xCol, yCol, *rows, zCol,
                                                  nBinsX, nBinsY, logFlag)
                    : QtFITSviewer::QtDensityPlot(customPlot, tab, xCol, yCol, zCol,
                                                  nBinsX, nBinsY, logFlag);
            }
        } catch (const std::exception& e) {
            QMessageBox::critical(nullptr, "Density error", QString::fromStdString(e.what()));
            return;
        }
        if (!colorMap) {
            QMessageBox::critical(nullptr, "Density error", "Density plot returned null.");
            return;
        }

        // A single complex or 2-vector column gives both axes, whichever selector holds it.
        rememberPlot(PlotSource::Kind::Density,
                     (*plottableY)[static_cast<size_t>((xIsNone ? yDensIdx : xDensIdx) - 1)].name,
                     (xIsNone || yIsNone) ? std::string() : (*plottableY)[static_cast<size_t>(yDensIdx - 1)].name,
                     rows);
        plotted->nBinsX = nBinsX;
        plotted->nBinsY = nBinsY;
        plotted->logX   = densLogX->isChecked();
        plotted->logY   = densLogY->isChecked();

        // QtDensityPlot only knows generic "x" and "y" labels, and has already drawn them.
        QString xLabel, yLabel;
        if (xIsNone || yIsNone) {
            const PlotCol& single = (*plottableY)[static_cast<size_t>((xIsNone ? yDensIdx : xDensIdx) - 1)];
            const DSL::dtype type = tab->getColumn(single.name)->getType();
            const bool isComplex = (type == DSL::dtype::tcplx || type == DSL::dtype::tdbcplx);
            const QString name = QString::fromStdString(single.name);
            const QString unit = single.unit.empty() ? QString() : " [" + QString::fromStdString(single.unit) + "]";
            xLabel = (isComplex ? "Re(" + name + ")" : name + "[0]") + unit;
            yLabel = (isComplex ? "Im(" + name + ")" : name + "[1]") + unit;
        } else {
            xLabel = colLabel((*plottableY)[static_cast<size_t>(xDensIdx - 1)]);
            yLabel = colLabel((*plottableY)[static_cast<size_t>(yDensIdx - 1)]);
        }
        // With the Log options the axes themselves hold log10 of the values.
        customPlot->xAxis->setLabel(densLogX->isChecked() ? "log10(" + xLabel + ")" : xLabel);
        customPlot->yAxis->setLabel(densLogY->isChecked() ? "log10(" + yLabel + ")" : yLabel);

        // With a weight, the cells hold that column's values rather than a count of rows.
        // A single complex or 2-vector column takes no weight, so its scale stays a count.
        const int weightIdx = zDensCombo->currentIndex();
        if (weightIdx > 0 && !xIsNone && !yIsNone && customPlot->plotLayout()->hasElement(0, 1)) {
            if (auto* colorScale = dynamic_cast<QCPColorScale*>(customPlot->plotLayout()->element(0, 1))) {
                const QString zLabel = colLabel((*plottableX)[static_cast<size_t>(weightIdx - 1)]);
                colorScale->axis()->setLabel(densLogZ->isChecked() ? "log10(" + zLabel + ")" : zLabel);
            }
        }

        if (!retain) {
            // No-retain path: keep the filled color map as-is, redrawn for the labels.
            customPlot->replot(QCustomPlot::rpQueuedReplot);
            markPlotted(rows);
            return;
        }

        // Retain path: extract the binned grid, remove the temporary color map,
        // and overlay contour lines with a fresh uniform color.

        auto* mapData = colorMap->data();
        const int gridW = mapData->keySize();
        const int gridH = mapData->valueSize();
        QCPRange xRange, yRange;
        xRange = mapData->keyRange();
        yRange = mapData->valueRange();

        // Build row-major Z grid: Z[row][col] = mapData->cell(col, row)
        std::vector<std::vector<double>> Z(
            static_cast<size_t>(gridH),
            std::vector<double>(static_cast<size_t>(gridW)));
        for (int row = 0; row < gridH; ++row)
            for (int col = 0; col < gridW; ++col)
                Z[static_cast<size_t>(row)][static_cast<size_t>(col)] =
                    mapData->cell(col, row);

        // Remove temporary density map and its color scale from the layout
        customPlot->removePlottable(colorMap.get());
        if (customPlot->plotLayout()->hasElement(0, 1)) {
            customPlot->plotLayout()->remove(customPlot->plotLayout()->element(0, 1));
            customPlot->plotLayout()->simplify();
        }

        // Compute contours (8 levels)
        QCustomContours contour(Z, 8);
        if (contour.lines().isEmpty()) {
            QMessageBox::warning(nullptr, "Density / Contour",
                "No contour levels could be computed (all bins may be equal or empty).");
            return;
        }

        // Each successive Retain+Plot call uses the next color in the palette
        const QColor cColor = contourPalette[(*contourColorIdx)++ % nContourColors];
        contour.setPen(QPen(cColor, 1.5));
        contour.clearColorGradient(); // uniform line color per spec
        contour.setLabelFontSize(0);  // suppress labels — clutters overlay

        // Render contour lines with coordinate remapping from grid-index space
        // to data space.  QCustomContours::render() uses raw grid indices [0, nBins)
        // as coordinates; we bypass it and remap manually so the contours align
        // with any scatter / histogram overlay already on the canvas.
        const double xScale = (gridW > 1)
            ? (xRange.upper - xRange.lower) / static_cast<double>(gridW - 1) : 1.0;
        const double yScale = (gridH > 1)
            ? (yRange.upper - yRange.lower) / static_cast<double>(gridH - 1) : 1.0;
        const double nanVal = std::numeric_limits<double>::quiet_NaN();

        for (const auto& lvlLines : contour.lines()) {
            const QVector<double>& xs = lvlLines.first;
            const QVector<double>& ys = lvlLines.second;
            auto* curve = new QCPCurve(customPlot->xAxis, customPlot->yAxis);
            curve->setPen(QPen(cColor, 1.5));
            curve->setScatterStyle(QCPScatterStyle::ssNone);
            for (int i = 0; i < xs.size(); ++i) {
                if (std::isnan(xs[i])) {
                    curve->addData(nanVal, nanVal);
                } else {
                    curve->addData(
                        xRange.lower + xs[i] * xScale,
                        yRange.lower + ys[i] * yScale);
                }
            }
        }

        customPlot->rescaleAxes();
        customPlot->replot(QCustomPlot::rpQueuedReplot);
        markPlotted(rows);
    };

    // --- doVectorRow lambda ---
    // Draws one row of a vector column. With a filter, a row it excludes moves to the first it keeps.
    auto doVectorRow = [=](const RowSet* rows) {
        if (vectorYCols->empty())
            return;

        size_t row = static_cast<size_t>(vecRowSpin->value() - 1);
        if (rows && !std::binary_search(rows->indices().begin(), rows->indices().end(), row)) {
            row = rows->indices().front();
            const QSignalBlocker block(vecRowSpin);
            vecRowSpin->setValue(static_cast<int>(row + 1));
        }

        const RowSet oneRow(std::vector<size_t>{row});
        const PlotCol& yColumn = (*vectorYCols)[static_cast<size_t>(vecYCombo->currentIndex())];
        const int xIdx = vecXCombo->currentIndex();

        std::vector<std::shared_ptr<QCPGraph>> graphs;
        try {
            graphs = (xIdx == 0)
                ? QtFITSviewer::QtScatter(customPlot, tab, yColumn.name, oneRow)
                : QtFITSviewer::QtScatter(customPlot, tab, (*vectorXCols)[static_cast<size_t>(xIdx - 1)].name, yColumn.name, oneRow);
        } catch (const std::exception& e) {
            QMessageBox::warning(nullptr, "Vector row", QString::fromStdString(e.what()));
            return;
        }

        const QColor col = scatterPalette[(*scatterColorIdx)++ % nScatterColors];
        for (auto& g : graphs) {
            g->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 4));
            g->setLineStyle(vecLineBox->isChecked() ? QCPGraph::lsLine : QCPGraph::lsNone);
            g->setPen(QPen(col));
        }

        customPlot->xAxis->setLabel((xIdx == 0) ? QString("Element") : colLabel((*vectorXCols)[static_cast<size_t>(xIdx - 1)]));
        customPlot->yAxis->setLabel(colLabel(yColumn));
        customPlot->rescaleAxes();
        customPlot->replot(QCustomPlot::rpQueuedReplot);
        *hasContent = true;
        showRowCount(1);
        showCategoryKey({}, {});
    };

    // --- Right panel: tab widget + Plot/Clear/Retain below it ---
    auto* rightPanel  = new QWidget;
    auto* rightLayout = new QVBoxLayout(rightPanel);
    rightLayout->setContentsMargins(4, 4, 4, 4);
    rightLayout->addWidget(tabWidget);
    rightLayout->addWidget(categoryBox);

    auto* plotBtn  = new QPushButton("Plot");
    auto* clearBtn = new QPushButton("Clear");
    auto* retainBox = new QCheckBox("Retain plot");
    auto* selectBtn = new QPushButton("Select region");
    selectBtn->setCheckable(true);
    selectBtn->setToolTip("While on, dragging on a scatter plot or a heatmap draws a box and lists "
                          "the rows inside it. Off, dragging moves the plot.");

    auto* btnLayout = new QHBoxLayout;
    btnLayout->addWidget(plotBtn);
    btnLayout->addWidget(clearBtn);
    rightLayout->addLayout(btnLayout);
    rightLayout->addWidget(selectBtn);
    rightLayout->addWidget(retainBox);
    rightPanel->setMinimumWidth(160);

    // --- doPlot lambda ---
    // Evaluates the row filter, then dispatches to the active tab's draw function.
    auto doPlot = [=]() {
        const std::string expression = filterEdit->text().trimmed().toStdString();

        // Before touching the canvas: a mistyped filter must leave the current plot in place.
        RowSet filtered;
        const RowSet* rows = nullptr;
        if (!expression.empty()) {
            try {
                filtered = tab->where(expression);
            } catch (const std::exception& e) {
                QMessageBox::warning(nullptr, "Row filter", QString::fromStdString(e.what()));
                return;
            }
            rows = &filtered;
        }

        // Same here: a bad epsilon or too many categories must leave the current plot in place.
        const int tabIdx = tabWidget->currentIndex();
        const int catIdx = categoryCombo->currentIndex();
        std::vector<TableCategory> categories;
        if (catIdx > 0 && (tabIdx == 0 || tabIdx == 2)) { // Scatter and Histogram only
            try {
                const auto& catCol = tab->getColumn((*categoryCols)[static_cast<size_t>(catIdx - 1)].name);

                double epsilon = 0.0;
                if (isRealColumn(*catCol)) {
                    const std::string text = epsilonEdit->text().trimmed().toStdString();
                    size_t consumed = 0;
                    epsilon = parseDouble(text, &consumed);
                    if (consumed != text.size() || !std::isfinite(epsilon) || epsilon < 0.0)
                        throw std::invalid_argument("Epsilon must be a positive number or 0, not '" + text + "'.");
                }

                categories = buildCategories(*catCol, rows ? *rows : tab->where(""), epsilon);
            } catch (const std::exception& e) {
                QMessageBox::warning(nullptr, "Category", QString::fromStdString(e.what()));
                return;
            }
        }

        const bool retain = retainBox->isChecked();
        highlightCategory(-1); // a retained plot must not stay dimmed
        if (!retain) clearCanvas();

        // A selection belongs to the plot it was made on.
        clearSelection();
        *plotted = PlotSource();

        // Nothing to draw, and the plotting kernels refuse empty data.
        if (rows && rows->empty()) {
            showRowCount(0);
            return;
        }

        switch (tabIdx) {
            case 0: doScatter(rows, categories);   break;
            case 1: doDensity(retain, rows);       break;
            case 2: doHistogram(rows, categories); break;
            case 3: doVectorRow(rows);             break;
            default: break;
        }
    };

    // --- presentSelection lambda ---
    // Lists the rows found on the plot, and rings the points they gave on a scatter plot.
    auto presentSelection = [=](const RowSet& rows, const PlotPoints& marks, const QString& description) {
        clearSelection();
        selectionLabel->setText(QString("%1 rows selected").arg(rows.size()));

        if (!marks.empty()) {
            QVector<double> xs, ys;
            xs.reserve(static_cast<int>(marks.size()));
            ys.reserve(static_cast<int>(marks.size()));
            for (const auto& mark : marks) {
                xs.push_back(mark.first);
                ys.push_back(mark.second);
            }
            *selectionMarks = customPlot->addGraph(customPlot->xAxis, customPlot->yAxis);
            (*selectionMarks)->removeFromLegend();
            (*selectionMarks)->setLayer("selectionMarks");
            (*selectionMarks)->setLineStyle(QCPGraph::lsNone);
            (*selectionMarks)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(QColor(255, 140, 0), 1.5), Qt::NoBrush, 10));
            (*selectionMarks)->setData(xs, ys);
        }
        customPlot->replot();

        if (!rows.empty())
            showSelectedRows(win, tab, rows, description);
    };

    // A box drawn on the plot: the rows with a point inside it.
    auto selectBox = [=](const QRect& box) {
        if (plotted->kind == PlotSource::Kind::None)
            return;

        const QRect r = box.normalized();
        const QCPRange xRange(std::min(customPlot->xAxis->pixelToCoord(r.left()), customPlot->xAxis->pixelToCoord(r.right())),
                              std::max(customPlot->xAxis->pixelToCoord(r.left()), customPlot->xAxis->pixelToCoord(r.right())));
        const QCPRange yRange(std::min(customPlot->yAxis->pixelToCoord(r.top()), customPlot->yAxis->pixelToCoord(r.bottom())),
                              std::max(customPlot->yAxis->pixelToCoord(r.top()), customPlot->yAxis->pixelToCoord(r.bottom())));
        try {
            PlotPoints marks;
            const RowSet rows = rowsInBox(*tab, *plotted, xRange, yRange, marks);
            // A heatmap keeps its cells visible: only a scatter plot gets rings.
            if (plotted->kind == PlotSource::Kind::Density)
                marks.clear();
            presentSelection(rows, marks, QString("box x %1 to %2, y %3 to %4")
                .arg(xRange.lower, 0, 'g', 6).arg(xRange.upper, 0, 'g', 6)
                .arg(yRange.lower, 0, 'g', 6).arg(yRange.upper, 0, 'g', 6));
        } catch (const std::exception& e) {
            QMessageBox::warning(nullptr, "Selection", QString::fromStdString(e.what()));
        }
    };

    // A click on the plot: the nearest scatter point, or the rows of a heatmap cell.
    auto selectAt = [=](const QPoint& position) {
        try {
            if (plotted->kind == PlotSource::Kind::Scatter) {
                PlotPoints marks;
                const RowSet rows = rowNearest(*tab, *plotted, customPlot->xAxis, customPlot->yAxis, position, 10.0, marks);
                if (!rows.empty()) // Clicking empty space does nothing.
                    presentSelection(rows, marks, QString("row %1").arg(rows.indices().front() + 1));
            } else if (plotted->kind == PlotSource::Kind::Density) {
                const double x = customPlot->xAxis->pixelToCoord(position.x());
                const double y = customPlot->yAxis->pixelToCoord(position.y());
                presentSelection(rowsInDensityCell(*tab, *plotted, x, y), PlotPoints(),
                                 QString("heatmap cell at x %1, y %2").arg(x, 0, 'g', 6).arg(y, 0, 'g', 6));
            }
        } catch (const std::exception& e) {
            QMessageBox::warning(nullptr, "Selection", QString::fromStdString(e.what()));
        }
    };

    // Previous / Next on the Vector row tab: through the filtered rows when a filter is set.
    auto stepVectorRow = [=](bool forward) {
        const std::string expression = filterEdit->text().trimmed().toStdString();
        const size_t current = static_cast<size_t>(vecRowSpin->value() - 1);

        if (expression.empty()) {
            vecRowSpin->setValue(vecRowSpin->value() + (forward ? 1 : -1));
            return;
        }

        RowSet rows;
        try {
            rows = tab->where(expression);
        } catch (const std::exception& e) {
            QMessageBox::warning(nullptr, "Row filter", QString::fromStdString(e.what()));
            return;
        }

        const auto& indices = rows.indices();
        if (forward) {
            const auto next = std::upper_bound(indices.begin(), indices.end(), current);
            if (next != indices.end())
                vecRowSpin->setValue(static_cast<int>(*next + 1));
        } else {
            const auto next = std::lower_bound(indices.begin(), indices.end(), current);
            if (next != indices.begin())
                vecRowSpin->setValue(static_cast<int>(*std::prev(next) + 1));
        }
    };

    // --- Central container: splitter only ---
    auto* splitter = new QSplitter(Qt::Horizontal);
    splitter->addWidget(leftPanel);
    splitter->addWidget(centerPanel);
    splitter->addWidget(rightPanel);
    splitter->setSizes({100, 600, 170});

    win->setCentralWidget(splitter);
    win->resize(1100, 680);

    // --- Wire connections ---
    // Plot button, or Return in the filter field
    QObject::connect(plotBtn, &QPushButton::clicked, win, [=]() { doPlot(); });
    QObject::connect(filterEdit, &QLineEdit::returnPressed, win, [=]() { doPlot(); });
    QObject::connect(epsilonEdit, &QLineEdit::returnPressed, win, [=]() { doPlot(); });

    // Categories only apply to Scatter and Histogram, a selection box to Scatter and Density
    QObject::connect(tabWidget, &QTabWidget::currentChanged, win, [=](int index) {
        categoryBox->setEnabled(index == 0 || index == 2);
        selectBtn->setEnabled(index == 0 || index == 1);
        if (!selectBtn->isEnabled())
            selectBtn->setChecked(false);
    });

    // Select region: dragging draws a box instead of moving the plot
    QObject::connect(selectBtn, &QPushButton::toggled, win, [=](bool on) {
        customPlot->setSelectionRectMode(on ? QCP::srmCustom : QCP::srmNone);
    });
    QObject::connect(customPlot->selectionRect(), &QCPSelectionRect::accepted, win, [=](const QRect& box, QMouseEvent*) {
        selectBox(box);
    });

    // A click is a press and release at the same place; a drag moves the plot or draws a box.
    auto pressedAt = std::make_shared<QPoint>();
    QObject::connect(customPlot, &QCustomPlot::mousePress, win, [=](QMouseEvent*) {
        *pressedAt = customPlot->mapFromGlobal(QCursor::pos());
    });
    QObject::connect(customPlot, &QCustomPlot::mouseRelease, win, [=](QMouseEvent* event) {
        const QPoint position = customPlot->mapFromGlobal(QCursor::pos());
        if (event->button() == Qt::LeftButton &&
            (position - *pressedAt).manhattanLength() <= 3 &&
            customPlot->axisRect()->rect().contains(position))
            selectAt(position);
    });

    // Vector row: a change of row redraws at once
    QObject::connect(vecRowSpin, QOverload<int>::of(&QSpinBox::valueChanged), win, [=](int) {
        if (tabWidget->currentIndex() == 3)
            doPlot();
    });
    QObject::connect(vecPrevBtn, &QPushButton::clicked, win, [=]() { stepVectorRow(false); });
    QObject::connect(vecNextBtn, &QPushButton::clicked, win, [=]() { stepVectorRow(true); });

    // Epsilon only applies to float and double columns
    QObject::connect(categoryCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), win, [=](int index) {
        epsilonEdit->setEnabled(index > 0
            && isRealColumn(*tab->getColumn((*categoryCols)[static_cast<size_t>(index - 1)].name)));
    });

    // Double click on a column of the left panel: window with its values and statistics
    QObject::connect(colTable, &QTableWidget::cellDoubleClicked, win, [=](int row, int) {
        const FITStable::clist columns = tab->listColumns();
        const auto& c = columns[static_cast<size_t>(row)];
        showColumnWindow(win, tab, c[0], c[1]);
    });

    // New column: computed in the table, then offered everywhere a column of its type can be chosen.
    // The filter and the double click read the table itself. Appending to the combos leaves
    // their current selections, and the index each one maps to in its list, where they are.
    auto addComputedColumn = [=]() {
        const QString definition = newColumnEdit->text().trimmed();
        if (definition.isEmpty())
            return;

        const FITSform* computed = nullptr;
        try {
            computed = &table->computeColumn(definition.toStdString());
        } catch (const std::exception& e) {
            QMessageBox::warning(nullptr, "New column", QString::fromStdString(e.what()));
            return;
        }

        const auto& column = tab->getColumn(computed->getName());
        const PlotCol entry{computed->getName(), computed->getUnit()};
        const QString label = colLabel(entry);

        if (isPlottableAsX(column)) {
            plottableX->push_back(entry);
            xScatterCombo->addItem(label);
            zDensCombo->addItem(label);
            histCombo->addItem(label);
        }
        if (isPlottableAsY(column)) {
            plottableY->push_back(entry);
            yScatterCombo->addItem(label);
            xDensCombo->addItem(label);
            yDensCombo->addItem(label);
        }
        if (isCategoryColumn(*column)) {
            categoryCols->push_back(entry);
            categoryCombo->addItem(label);
        }
        if (computed->getNelem() > 1 && isPlottableAsY(column)) {
            vectorYCols->push_back(entry);
            vecYCombo->addItem(label);
            tabWidget->setTabEnabled(3, true);
        }
        if (computed->getNelem() > 1 && isPlottableAsX(column)) {
            vectorXCols->push_back(entry);
            vecXCombo->addItem(label);
        }

        // Same row order as listColumns(), which the double click relies on.
        const int row = colTable->rowCount();
        colTable->insertRow(row);
        auto* nameItem = new QTableWidgetItem(QString::fromStdString(entry.name));
        QFont nameFont = nameItem->font();
        nameFont.setItalic(true);
        nameItem->setFont(nameFont);
        nameItem->setToolTip("Computed, in memory only: " + definition);
        colTable->setItem(row, 0, nameItem);
        colTable->setItem(row, 1, new QTableWidgetItem(QString::fromStdString(computed->getTTYPE())));
        colTable->scrollToItem(nameItem);
        columnCountLabel->setText(QString("Columns: %1").arg(colTable->rowCount()));

        newColumnEdit->clear();
    };
    QObject::connect(newColumnBtn, &QPushButton::clicked, win, [=]() { addComputedColumn(); });
    QObject::connect(newColumnEdit, &QLineEdit::returnPressed, win, [=]() { addComputedColumn(); });

    // Clear button: unconditional canvas wipe
    QObject::connect(clearBtn, &QPushButton::clicked, win, [=]() {
        clearCanvas();
    });

    // No auto-draw on combo change: the canvas stays blank until Plot is clicked.

    addCommandBar(win, tab->HDU());
    return win;
}

// ===========================================================================
// Colour compositing (opt-in via the HDU-list window's "Color" menu)
// ===========================================================================

// One image HDU usable as a colour channel.
struct ImageHduItem { int index; int naxis; int naxis3; std::string label; };

// List the image HDUs (NAXIS>=2) of a FITS file, with a human-readable label.
static std::vector<ImageHduItem> imageHdusOf(const std::string& path)
{
    std::vector<ImageHduItem> out;
    try {
        DataFile fm(path);
        std::vector<HduInfo> hdus = enumerateHDUs(fm);
        fm.Close();
        for (const auto& h : hdus) {
            if (h.type == HduInfo::Type::IMG && h.naxis >= 2) {
                std::string lbl = "HDU " + std::to_string(h.index) + " — " + h.extname
                                + " (" + std::to_string(h.naxis) + "D)";
                out.push_back({h.index, h.naxis, h.naxis3, lbl});
            }
        }
    } catch (...) { /* unreadable file → empty list */ }
    return out;
}

// Load one channel as a 2-D FITScube. For a 3-D HDU, extract layer (>=0).
static std::shared_ptr<const FITScube>
resolveChannelCube(const std::string& file, int hduIndex, int layer)
{
    DataFile fm(file);
    std::shared_ptr<FITScube> img = fm.GetImageAtIndex(hduIndex);
    if (img->GetDimension() == 3 && layer >= 0)
        img = img->Layer(static_cast<size_t>(layer));
    return std::const_pointer_cast<const FITScube>(img);
}

// Gradient presets offered per channel (RGB mode only).
enum class GradPreset { DefaultHue = 0, Viridis = 1, Grayscale = 2 };

static QCPColorGradient gradientFor(GradPreset p, int ch)
{
    if (p == GradPreset::Viridis)
        return QtFITS::makeViridisGradient();

    QCPColorGradient g;
    g.clearColorStops();
    if (p == GradPreset::Grayscale) {
        g.setColorStops({{0.0, QColor(0,0,0)}, {1.0, QColor(255,255,255)}});
        return g;
    }
    // DefaultHue: match QCPColorMapRGB's per-channel default ramp.
    switch (ch) {
        case 0:  g.setColorStops({{0.0, QColor(0,0,0)}, {1.0, QColor(255,0,0)}}); break;
        case 1:  g.setColorStops({{0.0, QColor(0,0,0)}, {1.0, QColor(0,255,0)}}); break;
        case 2:  g.setColorStops({{0.0, QColor(0,0,0)}, {1.0, QColor(0,0,255)}}); break;
        default: g.setColorStops({{0.0, QColor(0,0,0)}, {1.0, QColor(255,255,255)}}); break;
    }
    return g;
}

// Resolved output of the compose dialog.
struct ComposeResult {
    std::vector<std::shared_ptr<const FITScube>> channels;
    std::vector<QCPColorGradient>                gradients;   // per channel (RGB)
    std::vector<std::string>                     labels;      // per channel, for the dock
    std::string          file;                                 // of channel 0, offered for catalogues
    bool                 cmyk      = false;
    QCPAxis::ScaleType   scaleType = QCPAxis::stLinear;        // initial, per channel
};

// ---------------------------------------------------------------------------
// ColorComposeDialog — choose a colour model and assign FITS data to channels.
// ---------------------------------------------------------------------------
class ColorComposeDialog : public QDialog
{
public:
    // files[0] is the "current" file (used by the single-cube panel and as the
    // default for channel 0). Additional files seed the per-channel file combos
    // so a multi-file launch can map one file per channel out of the box.
    explicit ColorComposeDialog(const std::vector<std::string>& files,
                                QWidget* parent)
        : QDialog(parent),
          m_currentFile(files.empty() ? std::string() : files.front())
    {
        setWindowTitle(tr("Compose color image"));
        for (const std::string& f : files)
            if (std::find(m_files.begin(), m_files.end(), f) == m_files.end())
                m_files.push_back(f);
        if (m_files.empty()) m_files.push_back(m_currentFile);

        auto* root = new QVBoxLayout(this);

        // --- top controls -------------------------------------------------
        auto* form = new QFormLayout;
        m_modelCombo = new QComboBox(this);
        m_modelCombo->addItem("RGB (3 channels)");
        m_modelCombo->addItem("CMYK (4 channels)");
        form->addRow(tr("Color model:"), m_modelCombo);

        m_sourceCombo = new QComboBox(this);
        m_sourceCombo->addItem("From HDUs / files");
        m_sourceCombo->addItem("From single 3-D cube");
        form->addRow(tr("Source:"), m_sourceCombo);

        m_stretchCombo = new QComboBox(this);
        m_stretchCombo->addItem("Linear");
        m_stretchCombo->addItem("Logarithmic");
        form->addRow(tr("Stretch:"), m_stretchCombo);
        root->addLayout(form);

        // --- stacked source panels ---------------------------------------
        m_stack = new QStackedWidget(this);
        m_stack->addWidget(buildChannelPanel());   // page 0
        m_stack->addWidget(buildCubePanel());       // page 1
        root->addWidget(m_stack);

        // --- buttons ------------------------------------------------------
        auto* buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
        root->addWidget(buttons);
        connect(buttons, &QDialogButtonBox::accepted, this, &ColorComposeDialog::onAccept);
        connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);

        connect(m_modelCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, [this](int){ onModelChanged(); });
        connect(m_sourceCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, [this](int idx){ m_stack->setCurrentIndex(idx); });

        onModelChanged();   // build the initial channel rows
        resize(640, 360);
    }

    const ComposeResult& result() const { return m_result; }

private:
    struct ChannelRow {
        QLabel*      label     = nullptr;
        QComboBox*   fileCombo = nullptr;
        QPushButton* browseBtn = nullptr;
        QComboBox*   hduCombo  = nullptr;
        QSpinBox*    layerSpin = nullptr;
        QComboBox*   gradCombo = nullptr;
    };

    // ----- channel-table panel ----------------------------------------
    QWidget* buildChannelPanel()
    {
        m_channelPanel = new QWidget(this);
        m_channelGrid  = new QGridLayout(m_channelPanel);
        m_channelGrid->addWidget(new QLabel(tr("Channel")), 0, 0);
        m_channelGrid->addWidget(new QLabel(tr("File")),    0, 1);
        // col 2 = Browse button (no header)
        m_channelGrid->addWidget(new QLabel(tr("HDU")),     0, 3);
        m_channelGrid->addWidget(new QLabel(tr("Layer")),   0, 4);
        m_channelGrid->addWidget(new QLabel(tr("Gradient")),0, 5);
        return m_channelPanel;
    }

    void buildChannelRows(int n, const QStringList& names)
    {
        for (ChannelRow& r : m_rows) {
            delete r.label; delete r.fileCombo; delete r.browseBtn;
            delete r.hduCombo; delete r.layerSpin; delete r.gradCombo;
        }
        m_rows.clear();

        for (int i = 0; i < n; ++i) {
            ChannelRow r;
            r.label     = new QLabel(names.value(i), m_channelPanel);
            r.fileCombo = new QComboBox(m_channelPanel);
            r.browseBtn = new QPushButton(tr("Browse…"), m_channelPanel);
            r.hduCombo  = new QComboBox(m_channelPanel);
            r.layerSpin = new QSpinBox(m_channelPanel);
            r.layerSpin->setMinimum(0);
            r.layerSpin->setEnabled(false);
            r.gradCombo = new QComboBox(m_channelPanel);
            r.gradCombo->addItem("Default hue");
            r.gradCombo->addItem("Viridis");
            r.gradCombo->addItem("Grayscale");

            const int row = i + 1;
            m_channelGrid->addWidget(r.label,     row, 0);
            m_channelGrid->addWidget(r.fileCombo, row, 1);
            m_channelGrid->addWidget(r.browseBtn, row, 2);
            m_channelGrid->addWidget(r.hduCombo,  row, 3);
            m_channelGrid->addWidget(r.layerSpin, row, 4);
            m_channelGrid->addWidget(r.gradCombo, row, 5);

            const int idx = static_cast<int>(m_rows.size());
            populateFileCombo(r.fileCombo);
            // Default mapping when several files were supplied: file i → channel i
            // (e.g. file1→Red, file2→Green, file3→Blue). Set before connecting
            // signals so it does not fire onFileChanged prematurely.
            if (i < static_cast<int>(m_files.size()))
                r.fileCombo->setCurrentIndex(i);
            // File browsing is driven by a dedicated button (a clean top-level click),
            // NOT from inside the combo's currentIndexChanged — launching a modal file
            // dialog while the combo popup is closing prevented it from appearing.
            connect(r.browseBtn, &QPushButton::clicked,
                    this, [this, idx]{ onBrowse(idx); });
            connect(r.fileCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
                    this, [this, idx](int){ onFileChanged(idx); });
            connect(r.hduCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
                    this, [this, idx](int){ onHduChanged(idx); });

            m_rows.push_back(r);
            populateHduCombo(idx);
        }
    }

    void populateFileCombo(QComboBox* combo)
    {
        const QString keep = combo->currentData().toString();
        QSignalBlocker block(combo);
        combo->clear();
        for (const std::string& f : m_files) {
            QString full = QString::fromStdString(f);
            QString base = full.section('/', -1);
            combo->addItem(base, full);
        }
        int found = combo->findData(keep);
        combo->setCurrentIndex(found >= 0 ? found : 0);
    }

    // Selecting a different already-known file just refreshes that row's HDU list.
    void onFileChanged(int rowIdx)
    {
        populateHduCombo(rowIdx);
    }

    // Dedicated "Browse…" button handler — a clean top-level user event, so the
    // modal file dialog reliably appears (unlike launching it from a combo signal).
    void onBrowse(int rowIdx)
    {
        ChannelRow& r = m_rows[static_cast<size_t>(rowIdx)];

        // Start in the current file's directory so the user is not dropped in cwd.
        const QString startDir =
            QFileInfo(QString::fromStdString(m_currentFile)).absolutePath();

        // Use Qt's OWN file dialog on every platform (DontUseNativeDialog).
        // This is not a degraded mode — it is Qt's complete file browser that
        // navigates the whole filesystem identically on macOS and Linux. The
        // platform "native" panels were unreliable here (macOS returned empty;
        // a missing Linux portal/theme plugin can do the same).
        QFileDialog fd(this, tr("Open FITS file"), startDir);
        fd.setOption(QFileDialog::DontUseNativeDialog, true);
        fd.setFileMode(QFileDialog::ExistingFile);
        fd.setAcceptMode(QFileDialog::AcceptOpen);
        fd.setNameFilters({
            tr("FITS files (*.fits *.fit *.fts *.fz *.fits.gz *.fit.gz *.fts.gz)"),
#ifdef __WITH_ASDF__
            tr("ASDF files (*.asdf)"),
#endif
            tr("All files (*)")
        });

        if (fd.exec() != QDialog::Accepted) return;
        const QStringList sel = fd.selectedFiles();
        if (sel.isEmpty()) return;
        const QString picked = sel.first();

        const std::string p = picked.toStdString();
        if (std::find(m_files.begin(), m_files.end(), p) == m_files.end())
            m_files.push_back(p);
        for (ChannelRow& rr : m_rows) populateFileCombo(rr.fileCombo);

        const int found = r.fileCombo->findData(picked);
        if (found >= 0) r.fileCombo->setCurrentIndex(found);  // triggers populateHduCombo
        else            populateHduCombo(rowIdx);
    }

    void populateHduCombo(int rowIdx)
    {
        ChannelRow& r = m_rows[static_cast<size_t>(rowIdx)];
        const std::string file = r.fileCombo->currentData().toString().toStdString();
        QSignalBlocker block(r.hduCombo);
        r.hduCombo->clear();
        if (file.empty()) return;
        for (const ImageHduItem& it : imageHdusOf(file)) {
            r.hduCombo->addItem(QString::fromStdString(it.label), it.index);
            r.hduCombo->setItemData(r.hduCombo->count() - 1, it.naxis3, Qt::UserRole + 1);
        }
        // default each RGB channel to a distinct HDU when possible
        if (r.hduCombo->count() > rowIdx)
            r.hduCombo->setCurrentIndex(rowIdx);
        onHduChanged(rowIdx);
    }

    void onHduChanged(int rowIdx)
    {
        ChannelRow& r = m_rows[static_cast<size_t>(rowIdx)];
        const int naxis3 = r.hduCombo->currentData(Qt::UserRole + 1).toInt();
        const bool is3D = naxis3 > 0;
        r.layerSpin->setEnabled(is3D);
        r.layerSpin->setMaximum(is3D ? naxis3 - 1 : 0);
    }

    // ----- single-cube panel ------------------------------------------
    QWidget* buildCubePanel()
    {
        auto* panel = new QWidget(this);
        auto* form  = new QFormLayout(panel);
        m_cubeCombo = new QComboBox(panel);
        for (const ImageHduItem& it : imageHdusOf(m_currentFile))
            if (it.naxis >= 3)
                m_cubeCombo->addItem(QString::fromStdString(it.label), it.index);
        form->addRow(tr("3-D cube HDU:"), m_cubeCombo);

        m_cubeLayer[0] = new QSpinBox(panel); m_cubeLayer[0]->setValue(0);
        m_cubeLayer[1] = new QSpinBox(panel); m_cubeLayer[1]->setValue(1);
        m_cubeLayer[2] = new QSpinBox(panel); m_cubeLayer[2]->setValue(2);
        form->addRow(tr("Red   ← layer:"), m_cubeLayer[0]);
        form->addRow(tr("Green ← layer:"), m_cubeLayer[1]);
        form->addRow(tr("Blue  ← layer:"), m_cubeLayer[2]);

        auto updateMax = [this]() {
            const int naxis3 = m_cubeCombo->count()
                ? imageHdusOfNaxis3(m_cubeCombo->currentData().toInt()) : 0;
            for (auto* s : m_cubeLayer) s->setMaximum(naxis3 > 0 ? naxis3 - 1 : 0);
        };
        connect(m_cubeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, [updateMax](int){ updateMax(); });
        updateMax();
        return panel;
    }

    int imageHdusOfNaxis3(int hduIndex)
    {
        for (const ImageHduItem& it : imageHdusOf(m_currentFile))
            if (it.index == hduIndex) return it.naxis3;
        return 0;
    }

    void onModelChanged()
    {
        const bool cmyk = (m_modelCombo->currentIndex() == 1);
        if (cmyk) {
            buildChannelRows(4, {"Cyan", "Magenta", "Yellow", "Key (black)"});
            // CMYK uses fixed conversion + the channel table only.
            m_sourceCombo->setCurrentIndex(0);
            m_sourceCombo->setEnabled(false);
            m_stack->setCurrentIndex(0);
            for (ChannelRow& r : m_rows) { r.gradCombo->setEnabled(false); }
        } else {
            buildChannelRows(3, {"Red", "Green", "Blue"});
            m_sourceCombo->setEnabled(true);
            for (ChannelRow& r : m_rows) { r.gradCombo->setEnabled(true); }
        }
    }

    void onAccept()
    {
        const bool cmyk = (m_modelCombo->currentIndex() == 1);
        const bool cubeMode = (!cmyk && m_sourceCombo->currentIndex() == 1);

        ComposeResult res;
        res.cmyk      = cmyk;
        res.scaleType = (m_stretchCombo->currentIndex() == 1)
                            ? QCPAxis::stLogarithmic : QCPAxis::stLinear;

        try {
            if (cubeMode) {
                if (m_cubeCombo->count() == 0)
                    throw std::runtime_error("No 3-D image HDU found in this file.");
                const int hduIndex = m_cubeCombo->currentData().toInt();
                DataFile fm(m_currentFile);
                std::shared_ptr<FITScube> cube = fm.GetImageAtIndex(hduIndex);
                if (cube->GetDimension() < 3)
                    throw std::runtime_error("Selected HDU is not a 3-D cube.");
                res.file = m_currentFile;
                for (int ch = 0; ch < 3; ++ch) {
                    auto layer = cube->Layer(static_cast<size_t>(m_cubeLayer[ch]->value()));
                    res.channels.push_back(std::const_pointer_cast<const FITScube>(layer));
                    res.gradients.push_back(gradientFor(GradPreset::DefaultHue, ch));
                    res.labels.push_back("layer " + std::to_string(m_cubeLayer[ch]->value()));
                }
            } else {
                for (size_t ch = 0; ch < m_rows.size(); ++ch) {
                    ChannelRow& r = m_rows[ch];
                    const std::string file = r.fileCombo->currentData().toString().toStdString();
                    if (file.empty() || r.hduCombo->count() == 0)
                        throw std::runtime_error("Channel " + std::to_string(ch + 1)
                                                 + ": no file/HDU selected.");
                    const int hduIndex = r.hduCombo->currentData().toInt();
                    const int layer = r.layerSpin->isEnabled() ? r.layerSpin->value() : -1;
                    res.channels.push_back(resolveChannelCube(file, hduIndex, layer));
                    if (ch == 0)
                        res.file = file;
                    res.gradients.push_back(
                        gradientFor(static_cast<GradPreset>(r.gradCombo->currentIndex()),
                                    static_cast<int>(ch)));
                    std::string lbl = r.fileCombo->currentText().toStdString()
                                    + ":" + std::to_string(hduIndex);
                    if (layer >= 0) lbl += "/L" + std::to_string(layer);
                    res.labels.push_back(lbl);
                }
            }

            // Validate identical pixel dimensions across all channels.
            const size_t nx = res.channels[0]->Size(1);
            const size_t ny = res.channels[0]->Size(2);
            for (const auto& c : res.channels)
                if (c->Size(1) != nx || c->Size(2) != ny)
                    throw std::runtime_error("All channels must share the same pixel dimensions.");
        }
        catch (const std::exception& e) {
            QMessageBox::warning(this, tr("Compose color image"), e.what());
            return;
        }

        m_result = std::move(res);
        accept();
    }

    std::string                 m_currentFile;
    std::vector<std::string>    m_files;
    QComboBox*                  m_modelCombo   = nullptr;
    QComboBox*                  m_sourceCombo  = nullptr;
    QComboBox*                  m_stretchCombo = nullptr;
    QStackedWidget*             m_stack        = nullptr;
    QWidget*                    m_channelPanel = nullptr;
    QGridLayout*                m_channelGrid  = nullptr;
    std::vector<ChannelRow>     m_rows;
    QComboBox*                  m_cubeCombo    = nullptr;
    QSpinBox*                   m_cubeLayer[3] = {nullptr, nullptr, nullptr};
    ComposeResult               m_result;
};

// QDoubleSpinBox that accepts and displays scientific / engineering notation
// (e.g. "1e-5", "2.3E8"), over an essentially unlimited magnitude range.
// A stock QDoubleSpinBox only handles fixed-decimal text and would round small
// magnitudes to zero; here validate()/valueFromText()/textFromValue() are
// overridden and the decimals are kept very high so setValue() never rounds.
class SciDoubleSpinBox : public QDoubleSpinBox
{
public:
    explicit SciDoubleSpinBox(QWidget* parent = nullptr) : QDoubleSpinBox(parent)
    {
        setLocale(QLocale::c());            // '.' decimal point, no group separators
        setRange(-1e300, 1e300);
        setDecimals(1000);                  // visible text comes from textFromValue()
        m_validator.setNotation(QDoubleValidator::ScientificNotation);
        m_validator.setLocale(QLocale::c());
        m_validator.setRange(-1e300, 1e300, 1000);
    }

protected:
    QValidator::State validate(QString& text, int& pos) const override
    {
        return m_validator.validate(text, pos);
    }

    double valueFromText(const QString& text) const override
    {
        return text.toDouble();
    }

    QString textFromValue(double value) const override
    {
        // 'g' → shortest of decimal / scientific with 6 significant digits.
        return QString::number(value, 'g', 6);
    }

private:
    QDoubleValidator m_validator;
};

// Build the interactive colour-image viewer window from a resolved ComposeResult.
static QWidget* makeColorViewer(const ComposeResult& res)
{
    auto* mw = new QtFITSviewer(400, 250, 1100, 800);
    mw->setWindowTitle(QString("Color composite (%1 channels%2)")
                        .arg(static_cast<int>(res.channels.size()))
                        .arg(res.cmyk ? ", CMYK" : ""));
    QCustomPlot* plot = mw->GetPlot();

    // Render with raw data (OPTION=""); stretch is driven live by scale type below.
    QCPColorMapRGB* rgb = QtFITSviewer::OverlayChannels(plot->axisRect(), res.channels, res.cmyk);
    if (!rgb) { delete mw; return nullptr; }

    for (size_t ch = 0; ch < res.channels.size(); ++ch) {
        rgb->setChannelScaleType(static_cast<int>(ch), res.scaleType);
        if (!res.cmyk && ch < res.gradients.size())
            rgb->setChannelGradient(static_cast<int>(ch), res.gradients[ch]);
    }
    rgb->invalidate();

    plot->xAxis->setLabel("pixel [px]");
    plot->yAxis->setLabel("pixel [px]");
    plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    QCPAxisRect* image = plot->axisRect(0);
    image->axis(QCPAxis::atLeft)->setScaleRatio(image->axis(QCPAxis::atBottom), 1.0);
    QObject::connect(plot, &QCustomPlot::beforeReplot, [image]() {
        image->axis(QCPAxis::atLeft)->setScaleRatio(image->axis(QCPAxis::atBottom), 1.0);
    });

    // --- per-channel colour bars + control dock --------------------------
    // Each channel gets its own QCPColorScale beside the image (a dynamic,
    // draggable colour bar, like the monochrome viewer) AND a dock group with
    // typed min/max, stretch and gradient.  The colour bar is the single source
    // of truth for the range: typing in the spin boxes sets the bar's data
    // range, and the bar's dataRangeChanged signal pushes the range into the
    // QCPColorMapRGB (and back into the spin boxes).
    QCPMarginGroup* barMargins = new QCPMarginGroup(plot);
    image->setMarginGroup(QCP::msTop | QCP::msBottom, barMargins);

    // Colour ramp shown on each bar: the channel gradient in RGB mode, or a
    // single-hue C/M/Y/K reference ramp in CMYK mode.
    auto barGradient = [&res](int ch) -> QCPColorGradient {
        if (!res.cmyk && ch < static_cast<int>(res.gradients.size()))
            return res.gradients[ch];
        QCPColorGradient g; g.clearColorStops();
        if (ch >= 3) {
            // Key/black plate is subtractive: a higher value lays down MORE black
            // ink. Reverse the ramp so the bar reads white (low) → black (high),
            // matching the darkening effect (compositing itself is unchanged).
            g.setColorStops({{0.0, QColor(255, 255, 255)}, {1.0, QColor(0, 0, 0)}});
            return g;
        }
        QColor hue;
        switch (ch) {
            case 0:  hue = QColor(0, 255, 255); break;  // cyan
            case 1:  hue = QColor(255, 0, 255); break;  // magenta
            case 2:  hue = QColor(255, 255, 0); break;  // yellow
            default: hue = QColor(255, 255, 255); break;
        }
        g.setColorStops({{0.0, QColor(0, 0, 0)}, {1.0, hue}});
        return g;
    };

    auto* dock = new QDockWidget(QObject::tr("Channels"), mw);
    dock->setAllowedAreas(Qt::RightDockWidgetArea | Qt::LeftDockWidgetArea);
    auto* dockBody = new QWidget(dock);
    auto* dockLay  = new QVBoxLayout(dockBody);

    std::vector<QCPColorScale*> scales;
    std::vector<double>         boundLo, boundHi;

    for (size_t ch = 0; ch < res.channels.size(); ++ch) {
        const int   chi  = static_cast<int>(ch);
        const auto& cube = res.channels[ch];
        double lo = cube->GetMinimum();
        double hi = cube->GetMaximum();
        if (hi <= lo) hi = lo + 1.0;
        boundLo.push_back(lo);
        boundHi.push_back(hi);

        const QString chName = QString::fromStdString(
            ch < res.labels.size() ? res.labels[ch] : ("channel " + std::to_string(ch)));

        // ---- colour bar in the plot layout (column 1 + ch) -------------
        auto* cs = new QCPColorScale(plot);
        cs->setType(QCPAxis::atRight);
        cs->setGradient(barGradient(chi));
        cs->setDataRange(QCPRange(lo, hi));
        // No axis label (the dock group already names the channel). Fixed width +
        // compact number format so changing digit counts on the tick labels never
        // reflows the image / plot layout.
        cs->setMinimumSize(QSize(78, 0));
        cs->setMaximumSize(QSize(78, QWIDGETSIZE_MAX));
        cs->axis()->setNumberFormat("gb");
        cs->axis()->setNumberPrecision(3);
        plot->plotLayout()->addElement(0, 1 + chi, cs);
        cs->setMarginGroup(QCP::msTop | QCP::msBottom, barMargins);
        scales.push_back(cs);

        rgb->setChannelRange(chi, QCPRange(lo, hi));   // keep image in sync with bar

        // ---- dock group: typed min/max, stretch, gradient --------------
        auto* box  = new QGroupBox(chName, dockBody);
        auto* form = new QFormLayout(box);

        // Scientific/engineering-notation fields (accept e.g. "1e-5").
        QDoubleSpinBox* minSpin = new SciDoubleSpinBox(box);
        QDoubleSpinBox* maxSpin = new SciDoubleSpinBox(box);
        for (QDoubleSpinBox* s : {minSpin, maxSpin}) {
            // Commit on Enter / focus-out, NOT on every keystroke. With keyboard
            // tracking on, each digit fired valueChanged → set the bar range →
            // wrote the value back into the field mid-typing, which reset the
            // editor and rejected out-of-range entries before you finished typing.
            s->setKeyboardTracking(false);
            // Constant width so the dock does not reflow as the digit count changes.
            s->setFixedWidth(120);
        }
        minSpin->setValue(lo);
        maxSpin->setValue(hi);

        auto* stretchBox = new QComboBox(box);
        stretchBox->addItem("Linear");
        stretchBox->addItem("Logarithmic");
        stretchBox->setCurrentIndex(res.scaleType == QCPAxis::stLogarithmic ? 1 : 0);

        form->addRow(QObject::tr("min:"), minSpin);
        form->addRow(QObject::tr("max:"), maxSpin);
        form->addRow(QObject::tr("stretch:"), stretchBox);

        QComboBox* gradBox = nullptr;
        if (!res.cmyk) {
            gradBox = new QComboBox(box);
            gradBox->addItem("Default hue");
            gradBox->addItem("Viridis");
            gradBox->addItem("Grayscale");
            form->addRow(QObject::tr("gradient:"), gradBox);
        }
        dockLay->addWidget(box);

        // ---- wiring ----------------------------------------------------
        // typed min/max → colour-bar data range
        auto pushSpins = [cs, minSpin, maxSpin]() {
            double a = minSpin->value(), b = maxSpin->value();
            if (b <= a) b = a + 1e-9;
            cs->setDataRange(QCPRange(a, b));
        };
        QObject::connect(minSpin, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                         plot, [pushSpins](double){ pushSpins(); });
        QObject::connect(maxSpin, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                         plot, [pushSpins](double){ pushSpins(); });

        // colour-bar range (drag / zoom / typed) → channel range + spin sync
        QObject::connect(cs, &QCPColorScale::dataRangeChanged, plot,
            [rgb, plot, chi, minSpin, maxSpin](const QCPRange& rng) {
                rgb->setChannelRange(chi, rng);
                rgb->invalidate();
                { QSignalBlocker b1(minSpin); QSignalBlocker b2(maxSpin);
                  minSpin->setValue(rng.lower); maxSpin->setValue(rng.upper); }
                plot->replot();
            });

        // stretch → channel scale type (+ bar axis scale type when range > 0)
        QObject::connect(stretchBox, QOverload<int>::of(&QComboBox::currentIndexChanged),
            plot, [rgb, plot, cs, chi](int idx) {
                const bool log = (idx == 1);
                rgb->setChannelScaleType(chi, log ? QCPAxis::stLogarithmic : QCPAxis::stLinear);
                cs->axis()->setScaleType((log && cs->dataRange().lower > 0)
                                         ? QCPAxis::stLogarithmic : QCPAxis::stLinear);
                rgb->invalidate();
                plot->replot();
            });

        // gradient → channel gradient + bar gradient
        if (gradBox)
            QObject::connect(gradBox, QOverload<int>::of(&QComboBox::currentIndexChanged),
                plot, [rgb, plot, cs, chi](int idx) {
                    const QCPColorGradient g = gradientFor(static_cast<GradPreset>(idx), chi);
                    rgb->setChannelGradient(chi, g);
                    cs->setGradient(g);
                    rgb->invalidate();
                    plot->replot();
                });

        if (res.scaleType == QCPAxis::stLogarithmic && lo > 0)
            cs->axis()->setScaleType(QCPAxis::stLogarithmic);
    }
    dockLay->addStretch(1);

    // Image dominates the width; the colour bars stay narrow.
    plot->plotLayout()->setColumnStretchFactor(0, 8);
    for (size_t ch = 0; ch < scales.size(); ++ch)
        plot->plotLayout()->setColumnStretchFactor(1 + static_cast<int>(ch), 0.2);

    auto* scroll = new QScrollArea(dock);
    scroll->setWidgetResizable(true);
    scroll->setWidget(dockBody);
    dock->setWidget(scroll);
    mw->addDockWidget(Qt::RightDockWidgetArea, dock);

    // --- command bar: reset levels (restore each bar to its data min/max) -
    auto* toolbar = new QToolBar("Options", mw);
    toolbar->setMovable(false);
    toolbar->setFloatable(false);
    auto* resetBtn = new QPushButton("Reset levels", toolbar);
    resetBtn->setFlat(true);
    QObject::connect(resetBtn, &QPushButton::clicked, mw,
        [scales, boundLo, boundHi, plot]() {
            for (size_t ch = 0; ch < scales.size(); ++ch)
                scales[ch]->setDataRange(QCPRange(boundLo[ch], boundHi[ch]));
            plot->replot();
        });
    toolbar->addWidget(resetBtn);
    mw->addToolBar(Qt::TopToolBarArea, toolbar);

    // Channels only have to share their pixel grid: catalogues use the WCS of the first one.
    addCatalogueTools(mw, plot, res.channels[0], res.file);
    plot->setNoAntialiasingOnDrag(true);

    // --- status-bar pixel readout ----------------------------------------
    // Connect a no-arg lambda to mouseMove (the QMouseEvent* arg is discarded);
    // read the cursor via QCursor::pos() exactly like the monochrome mouseHasMove.
    std::vector<std::shared_ptr<const FITScube>> chans = res.channels;
    QObject::connect(plot, &QCustomPlot::mouseMove, plot,
        [plot, image, chans, mw]() {
            const QPoint cursor = plot->mapFromGlobal(QCursor::pos());
            if (plot->axisRectAt(cursor) != image) return;
            const double kx = image->axis(QCPAxis::atBottom)->pixelToCoord(cursor.x());
            const double ky = image->axis(QCPAxis::atLeft)->pixelToCoord(cursor.y());
            const long ix = static_cast<long>(kx + 0.5);
            const long iy = static_cast<long>(ky + 0.5);
            if (ix < 0 || iy < 0
                || ix >= static_cast<long>(chans[0]->Size(1))
                || iy >= static_cast<long>(chans[0]->Size(2))) {
                if (mw->statusBar()) mw->statusBar()->clearMessage();
                return;
            }
            QString msg = QString("[%1 ; %2] ").arg(ix).arg(iy);
            for (size_t ch = 0; ch < chans.size(); ++ch) {
                double v = 0.0;
                try { v = chans[ch]->DoubleValueAtPixel(
                          {static_cast<size_t>(ix), static_cast<size_t>(iy)}); }
                catch (...) {}
                msg += QString("ch%1=%2  ").arg(static_cast<int>(ch)).arg(v);
            }
            if (mw->statusBar()) mw->statusBar()->showMessage(msg);
        });

    mw->SetPlot(plot);
    mw->show();
    return mw;
}

// Open the compose dialog seeded with one or more files; on Accept, build and
// show the colour viewer. Returns the viewer window (or nullptr if cancelled /
// nothing composed). Invoked both by HduListWindow::colorComposeRequested (one
// file) and by the multi-file command-line launch path.
static QWidget* openColorComposeDialog(const std::vector<std::string>& files,
                                       QWidget* parent)
{
    if (files.empty()) return nullptr;
    ColorComposeDialog dlg(files, parent);
    if (dlg.exec() != QDialog::Accepted) return nullptr;
    if (dlg.result().channels.empty()) return nullptr;
    return makeColorViewer(dlg.result());   // top-level window; stays alive after return
}

static QWidget* makeViewer(const HduInfo& hdu,
                            const std::string& fitsPath,
                            const UOption& opt)
{
    // Cancel and join every background thread from the previous 2D viewer before
    // doing anything else.  This must happen here — not inside make2DViewer —
    // because switching to a header-only or table HDU would otherwise leave the
    // old threads alive and calling invokeMethod on a deleted QCustomPlot.
    _cancel_wcs = true;
    if (_draw_histo.joinable()) _draw_histo.join();
    if (_draw_grid.joinable())  _draw_grid.join();
    if (_print_info.joinable()) _print_info.join();
    _cancel_wcs = false;

    if (hdu.type == HduInfo::Type::IMG) {
        switch (hdu.naxis) {
            case 0:
                return makeHeaderOnlyViewer(hdu, fitsPath);
            case 1:
                QMessageBox::warning(nullptr, "Not implemented",
                    QString("1D image viewer not yet implemented (HDU %1).").arg(hdu.index));
                return notYetImpl("1D image viewer", hdu);
            case 2:
                return make2DViewer(hdu, fitsPath, opt);
            case 3:
                return make3DViewer(hdu, fitsPath, opt);
            default:
                QMessageBox::warning(nullptr, "Unsupported HDU",
                    QString("NAXIS=%1 not supported (HDU %2).").arg(hdu.naxis).arg(hdu.index));
                return notYetImpl("N-D image viewer (NAXIS > 3)", hdu);
        }
    } else if (hdu.type == HduInfo::Type::TBL || hdu.type == HduInfo::Type::BTBL) {
        return makeFITSTableViewer(hdu, fitsPath);
    } else {
        QMessageBox::warning(nullptr, "Unsupported HDU",
            QString("Unknown HDU type (HDU %1).").arg(hdu.index));
        return notYetImpl("Unknown HDU type", hdu);
    }
}

#endif

int main(int argc, char **argv)
{
#ifdef __WITH_QT__

    try
    {
        UOption fOPT = UOption();
        if (!fOPT.InitParameters(argc, argv))
            return EXIT_FAILURE;

        std::ostream* fp = &std::cout;
        std::ofstream fout;
        if (fOPT.GetOutput().size() > 1)
        {
            fout.open(fOPT.GetOutput());
            fp = &fout;
        }

        int16_t debugLevel = fOPT.DebugLevel();
        if (debugLevel < 1)
            DSL::verbose = DSL::verboseLevel::VERBOSE_NONE;
        else if (debugLevel == 1)
            DSL::verbose = DSL::verboseLevel::VERBOSE_BASIC;
        else if (debugLevel == 2)
            DSL::verbose = DSL::verboseLevel::VERBOSE_DETAIL;
        else
            DSL::verbose = DSL::verboseLevel::VERBOSE_DEBUG;

        flist in_file = fOPT.GetFiles();

        if (in_file.size() < 1)
            throw std::runtime_error("Not enough arguments. Please provide a filename to open.");

        // Headless path: --draw disables GUI
        if (UOption::ddraw == UOption::drawing_mode::OFF)
        {
            if (in_file.size() > 1)
                throw std::runtime_error("Headless mode (--draw) handles a single file. "
                    "Provide one file, or drop --draw to compose a colour image from several.");

            DataFile ffile(in_file[0]);
            std::shared_ptr<FITScube> img = ffile.GetImageAtIndex(fOPT.Header());

            if (fOPT.ROI.size() > 1)
            {
                size_t X0 = 0, nX = img->Size(1);
                size_t Y0 = 0, nY = img->Size(2);
                if (fOPT.ROI.find('x') != fOPT.ROI.end()) {
                    X0 = fOPT.ROI.find('x')->second.first;
                    nX = fOPT.ROI.find('x')->second.second;
                }
                if (fOPT.ROI.find('y') != fOPT.ROI.end()) {
                    Y0 = fOPT.ROI.find('y')->second.first;
                    nY = fOPT.ROI.find('y')->second.second;
                }
                auto tmp = img->Window(X0, Y0, nX, nY);
                std::swap(img, tmp);
            }

            double sum=0, mean=0, median=0, rms=0, min=0, max=0, pp25=0, pp75=0;

            auto fut_sum    = std::async(std::launch::async, [img]{ return img->GetSum(); });
            auto fut_mean   = std::async(std::launch::async, [img]{ return img->GetMean(); });
            auto fut_median = std::async(std::launch::async, [img]{ return img->GetMedian(); });
            auto fut_rms    = std::async(std::launch::async, [img]{ return img->GetVariance(); });
            auto fut_min    = std::async(std::launch::async, [img]{ return img->GetMinimum(); });
            auto fut_max    = std::async(std::launch::async, [img]{ return img->GetMaximum(); });
            auto fut_pp25   = std::async(std::launch::async, [img]{ return img->Getpercentil(0.05); });
            auto fut_pp75   = std::async(std::launch::async, [img]{ return img->Getpercentil(0.95); });

            sum    = fut_sum.get();
            mean   = fut_mean.get();
            median = fut_median.get();
            rms    = fut_rms.get();
            min    = fut_min.get();
            max    = fut_max.get();
            pp25   = fut_pp25.get();
            pp75   = fut_pp75.get();

            if ((DSL::verbose & DSL::verboseLevel::VERBOSE_BASIC) == DSL::verboseLevel::VERBOSE_BASIC)
            {
                QtFITS::msg_m.lock();
                std::cout << "\033[01m=====\033[34m" << in_file[0] << "\033[0m" << std::endl
                          << "\033[31m     |- AVG :\033[0m" << mean   << " " << ((img->HDU().Exists("BUNIT"))?img->HDU().GetValueForKey("BUNIT"):"") << std::endl
                          << "\033[31m     |- MED :\033[0m" << median << " " << ((img->HDU().Exists("BUNIT"))?img->HDU().GetValueForKey("BUNIT"):"") << std::endl
                          << "\033[31m     |- 25% :\033[0m" << pp25   << " " << ((img->HDU().Exists("BUNIT"))?img->HDU().GetValueForKey("BUNIT"):"") << std::endl
                          << "\033[31m     |- 75% :\033[0m" << pp75   << " " << ((img->HDU().Exists("BUNIT"))?img->HDU().GetValueForKey("BUNIT"):"") << std::endl;
                QtFITS::msg_m.unlock();
            }

            if (fOPT.ArrayInfo())
            {
                QtFITS::msg_m.lock();
                std::cout << "\033[32m==== " << in_file[0] << "[" << fOPT.Header() << "]\033[0m" << std::endl;
                std::cout << " FILE |  NAME | MEAN | MEDIAN | Q1 | Q3 | RMS | MIN | MAX | SUM" << std::endl;
                std::cout << fOPT.GetFiles()[0] << " | " << img->GetName() << " | "
                          << mean << " | " << median << " | " << pp25 << " | " << pp75
                          << " | " << rms << " | " << min << " | " << max << " | " << sum << std::endl;
                QtFITS::msg_m.unlock();
            }

            if (fOPT.doPrint())
            {
                if (fOPT.GetOutput().size() > 0)
                    std::cout << " ... WRITING TO " << fOPT.GetOutput() << " plz, wait!" << std::endl;
                img->Print();
            }

            img.reset();
            ffile.Close();
            return EXIT_SUCCESS;
        }

        // GUI path
        QApplication a(argc, argv);

        // macOS native message boxes show the application's icon, a generic one for this
        // bare executable; Qt's own show the standard warning and error signs.
        QApplication::setAttribute(Qt::AA_DontUseNativeDialogs);

        // Several files on the command line → go straight to the colour compose
        // dialog (seeded with every file, mapped one-per-channel by default)
        // instead of treating it as an error.
        if (in_file.size() > 1)
        {
            std::vector<std::string> files(in_file.begin(), in_file.end());
            QWidget* viewer = openColorComposeDialog(files, nullptr);
            if (!viewer) return EXIT_SUCCESS;   // user cancelled the dialog
            const int ret = a.exec();
            if (_draw_histo.joinable()) _draw_histo.join();
            if (_draw_grid.joinable())  _draw_grid.join();
            if (_print_info.joinable()) _print_info.join();
            if (_print_.joinable())     _print_.join();
            return ret;
        }

        DataFile ffile(in_file[0]);
        std::vector<HduInfo> hdus = enumerateHDUs(ffile);
        ffile.Close();

        if (hdus.empty())
        {
            QMessageBox::critical(nullptr, "FITS error", "No readable HDUs found in file.");
            return EXIT_FAILURE;
        }

        const std::string fitsPath = in_file[0];
        auto factory = [&fOPT, fitsPath](const HduInfo& hdu) -> QWidget* {
            return makeViewer(hdu, fitsPath, fOPT);
        };

        HduListWindow listWin(in_file[0], hdus, factory);

        // Opt-in colour compositing: the "Color ▸ Compose color image…" menu emits
        // colorComposeRequested(); open the channel-selection dialog in response.
        QObject::connect(&listWin, &HduListWindow::colorComposeRequested,
                         [fitsPath, &listWin]() {
                             openColorComposeDialog({fitsPath}, &listWin);
                         });

        listWin.show();

        if (fOPT.HduExplicit())
            listWin.openViewer(fOPT.Header());

        const int ret = a.exec();

        // Join background threads started by Draw() to avoid std::terminate
        // on global std::thread destructor when threads are still joinable.
        if (_draw_histo.joinable()) _draw_histo.join();
        if (_draw_grid.joinable())  _draw_grid.join();
        if (_print_info.joinable()) _print_info.join();
        if (_print_.joinable())     _print_.join();

        return ret;
    }
    catch (std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        throw;
    }
    catch (...)
    {
        return EXIT_FAILURE;
    }
#else
    std::cerr << "\033[43;34m*** [FitsViewer]: Errors *** " << std::endl
              << "THIS APPLICATION REQUIRE Qt5+ TO BE INSTALLED" << std::endl
              << "PLZ, INSTALL Qt5 or Qt6 AND RETRY\033[0m" << std::endl;
    return EXIT_FAILURE;
#endif

    return EXIT_SUCCESS;
}