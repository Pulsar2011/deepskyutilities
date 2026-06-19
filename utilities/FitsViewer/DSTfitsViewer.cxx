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
#include <DSTdisplay/QCustomContours.h>
#include <DSTdisplay/QtFitsPloting.h>

#endif

#include <limits>

#include <DSTCore/OptionManager.h>
#include <DSTCore/DSS_version.h>
#include <DSTools/DST_version.h>
#include <DSTdisplay/DSD_version.h>
#include <DSTfits/FITSmanager.h>
#include <DSTfits/FITSimg.h>
#include <DSTfits/FITShdu.h>
#include <DSTmath/DSTmath.h>


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

typedef const std::vector<std::string> flist;
typedef std::map<char,std::pair<size_t,size_t> > roi;

size_t histo_bins=100;
double minZscale=std::numeric_limits<double>::quiet_NaN();
double maxZscale=std::numeric_limits<double>::quiet_NaN();

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
        LOG   =0x1,
        QUAD  =0x2,
        SQRT  =0x3,
        OFF   =0xF,
    };
    
private:
    int16_t HDU;
    size_t RBIN;
    std::string opt_string;
    rmode MODE;
    bool do_Print;
    bool a_info;
    bool hdu_explicit;

public:
    
    static drawing_mode ddraw;
    roi ROI;
    
    UOption():DOption(),HDU(1),RBIN(0),opt_string(),MODE(rmode::NONE), do_Print(false),a_info(false),hdu_explicit(false),ROI()
    {
        AddDescription("Open a FITS header block and display data.");
        AddOption("--hdu",DParam::fparameter::uint8,"Header block to be open. If not given, the primary header will be opened");
        AddOption("--opt/++opt",DParam::fparameter::string,"Define the color ramp used to display the image. Possible value are \033[34m'sqrtz'\033[0m, \033[34m'logz'\033[0m, \033[34m'quadz'\033[0m or \033[34m'linearz'\033[32m (default)\033[0m.");
        AddOption("--bins/-b",DParam::fparameter::uint16,"Number of bin to display the histograms of the images.");
        AddOption("--roi",{DParam::fparameter::xoffset,DParam::fparameter::xsize,DParam::fparameter::yoffset,DParam::fparameter::ysize},"Define a refgion of interest to be ploted. Other data outside this ROI are disregarded. The ROI shall be given in pixel not in WCS");
        AddOption("--print","Print pixel content to stdout.");
        AddOption("--info","Print statistical information on the array distribution.");
        AddOption("--draw","Disable drawing!");
        AddOption("--zscale/-z",{DParam::fparameter::min,DParam::fparameter::max},"Select the min and max value for the color scale.");
        AddParam("Fits_Filename","FITS file URL and name that will be opened for analyses.");    
        
        DOption::prg_name += "FitsViewer -- "+DST::gGIT::version();

        if(DSS::gGIT::version() < "v2.0.0")
            throw std::runtime_error("[ERROR] UOption::UOption: Incompatible DSS version, please update your DSS installation to v2.0.0 at minimum.");
        if(DSD::gGIT::version() < "v1.2.2")
            throw std::runtime_error("[ERROR] UOption::UOption: Incompatible DSD version, please update your DSD installation to v1.2.2 at minimum.");
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
        else if(!it->compare("++opt"))
        {
            opt_string.clear();
            if(std::next(it) != arg.cend()){it++;}else{return false;}
            opt_string = std::string(*it);
            
            if(opt_string.find("sqrtz") != std::string::npos)
                ddraw = drawing_mode::SQRT;
            else if(opt_string.find("logz") != std::string::npos)
                ddraw = drawing_mode::LOG;
            else if(opt_string.find("quadz") != std::string::npos)
                ddraw = drawing_mode::QUAD;
            else
                ddraw = drawing_mode::LINEAR;
            
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
            opt_string.clear();
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
    inline const std::string&  draw()  const {return opt_string;}
    inline bool doPrint() const {return do_Print;}
    inline bool ArrayInfo() const {return a_info;}
    inline bool HduExplicit() const {return hdu_explicit;}

};

UOption::drawing_mode UOption::ddraw = UOption::drawing_mode::LINEAR;

#ifdef __WITH_QT__

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
        if(in->getNumberOfWCS() > 0 && in->GetDimension() == 2)
        {
            worldCoords wc = in->WorldCoordinates(std::vector<size_t>({ix, iy}));
            std::vector< std::string > ctype = in->getWCS().CTYPE(0);
            std::vector< std::string > cunit = in->getWCS().CUNIT(0);

            std::string xstr = (ctype[0].size() > 0)?ctype[0]:"\u03b1"; 
            if(xstr.find('-') != std::string::npos) xstr.resize(xstr.find_first_of('-'));
            std::string ystr = (ctype[1].size() > 0)?ctype[1]:"\u03b4"; 
            if(ystr.find('-') != std::string::npos) ystr.resize(ystr.find_first_of('-'));

            std::string xUnit = (cunit[0].size() > 0)?" "+cunit[0]:""; if(xUnit == " deg" || xUnit == " DEG" ) xUnit = "\u00b0";
            std::string yUnit = (cunit[1].size() > 0)?" "+cunit[1]:""; if(yUnit == " deg" || yUnit == " DEG" ) yUnit = "\u00b0";

            wcs_msg = QString("%1=%2%4  %5=%6%7 ").arg(QString::fromStdString(xstr)).arg(wc[0],0,'f',3).arg(QString::fromStdString(xUnit)).arg(QString::fromStdString(ystr)).arg(wc[1],0,'f',3).arg(QString::fromStdString(yUnit));
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

    if(in->getNumberOfWCS() < 1)
        return;

    if( in->getWCS().CRPIX(0,1) == 0 && in->getWCS().CRPIX(0,2) == 0 &&
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
                    lbl->setFont(QFont("Monospace", 7));
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

QCustomPlot* Draw(const std::shared_ptr<const FITScube>& in, std::string opt)
{
    if(UOption::ddraw == UOption::drawing_mode::OFF)
        return NULL;
    
    QtFITSviewer *mw = new QtFITSviewer(400, 250, 1200, 900);
    mw->setWindowTitle(in->GetName().c_str());
    mw->linkResizeEvent(onResize);

    //size_t xsize = std::min(static_cast<size_t>(900),in->Size(1)+200);
    //size_t ysize = std::min(static_cast<size_t>(900),in->Size(2)+200);

    QCustomPlot *customPlot = mw->GetPlot();
    
    QtFITSviewer::QtDisplay(customPlot,in,opt);
    
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
    FITSmanager fm(fitsPath);
    const FITShdu fitshdu = *fm.GetHeaderAtIndex(hdu.index);

    auto* win = new QMainWindow;
    win->setWindowTitle(QString("HDU %1 — %2 (header only)")
        .arg(hdu.index)
        .arg(QString::fromStdString(hdu.extname)));

    win->setCentralWidget(makeFitsHeaderTable(fitshdu, win));
    win->resize(700, 400);
    return win;
}

static QWidget* make2DViewer(const HduInfo& hdu,
                              const std::string& fitsPath,
                              const UOption& opt)
{
    // Thread cleanup is now handled by makeViewer before dispatching here.

    FITSmanager fm(fitsPath);
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
    QCustomPlot* plt = Draw(c_img, opt.draw());
    if (!plt) return nullptr;

    QMainWindow* mainWin = qobject_cast<QMainWindow*>(plt->window());
    if (mainWin)
        addCommandBar(mainWin, img->HDU(), plt);
    return mainWin ? mainWin : dynamic_cast<QWidget*>(plt->window());
}

static QWidget* make3DViewer(const HduInfo& hdu,
                              const std::string& fitsPath,
                              const UOption& opt)
{
    FITSmanager fm(fitsPath);
    std::shared_ptr<FITScube> cube = fm.GetImageAtIndex(hdu.index);
    auto c_cube = std::const_pointer_cast<const FITScube>(cube);
    int nLayers = hdu.naxis3;
    std::string drawOpt = opt.draw();

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

    QtFITSviewer::QtDisplayLayer(customPlot, c_cube, 1, drawOpt);

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
    auto* arrowStyle = new ArrowSliderStyle(slider->style());
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
                     [customPlot, c_cube, px_pdf, px_cdf, layerLabel, nLayers, nx, ny, drawOpt](int layer0) {
                         int layer1 = layer0 + 1;

                         // Preserve current color scale range and axis limits across layer changes
                         QCPRange savedRange, savedXRange, savedYRange;
                         bool hasRange = false;
                         if (customPlot->plotLayout()->hasElement(0, 1)) {
                             auto* cs = dynamic_cast<QCPColorScale*>(
                                 customPlot->plotLayout()->element(0, 1));
                             if (cs) { savedRange = cs->dataRange(); hasRange = true; }
                         }
                         QCPAxisRect* imageRect = customPlot->axisRect(0);
                         savedXRange = imageRect->axis(QCPAxis::atBottom)->range();
                         savedYRange = imageRect->axis(QCPAxis::atLeft)->range();

                         QtFITSviewer::QtDisplayLayer(customPlot, c_cube, static_cast<size_t>(layer1), drawOpt);

                         // Restore color scale range
                         if (hasRange) {
                             if (customPlot->plotLayout()->hasElement(0, 1)) {
                                 auto* cs = dynamic_cast<QCPColorScale*>(
                                     customPlot->plotLayout()->element(0, 1));
                                 if (cs) cs->setDataRange(savedRange);
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

static QWidget* makeFITSTableViewer(const HduInfo& hdu, const std::string& fitsPath)
{
    // --- Load the table ---
    FITSmanager fm(fitsPath);
    auto tab = std::const_pointer_cast<const FITStable>(fm.GetTableAtIndex(hdu.index));

    // --- Column metadata ---
    const FITStable::clist cols = tab->listColumns(); // each entry: {name, TFORM, unit}
    const size_t nCols = tab->ncols();
    const size_t nRows = tab->nrows();

    // Build separate plottable lists for X and Y combo boxes.
    // plottableX: real/integer scalar columns (no complex) — used for Scatter X, Density Z weight, Histogram.
    // plottableY: all numeric columns including complex — used for Scatter Y, Density X/Y.
    // Complex columns are allowed in Density X/Y selectors; when one axis is left at "-- (none) --"
    // and the other is complex, QtDensityPlot's single-column overload splits real/imag parts.
    struct PlotCol { std::string name; std::string unit; };
    std::vector<PlotCol> plottableX;
    std::vector<PlotCol> plottableY;

    for (const auto& entry : cols) {
        const std::string& name = entry[0];
        const std::string& unit = entry[2];
        const auto& col = tab->getColumn(name);
        if (isPlottableAsX(col)) plottableX.push_back({name, unit});
        if (isPlottableAsY(col)) plottableY.push_back({name, unit});
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
    leftLayout->addWidget(new QLabel(QString("Columns: %1").arg(nCols)));
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
    leftLayout->addWidget(colTable, 1);
    leftPanel->setMinimumWidth(120);

    // --- Centre panel: plot area ---
    auto* customPlot = new QCustomPlot;
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes);
    customPlot->axisRect()->setupFullAxesBox(true);
    customPlot->legend->setVisible(false);  // axis labels identify columns; legend causes layout crash on first show
    customPlot->setMinimumWidth(200);

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
    for (const auto& c : plottableX)
        xScatterCombo->addItem(colLabel(c));

    // Y: all numeric columns including complex; no sentinel — column required
    for (const auto& c : plottableY)
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
    for (const auto& c : plottableY) {
        xDensCombo->addItem(colLabel(c));
        yDensCombo->addItem(colLabel(c));
    }

    // Z/weight: "-- (none) --" sentinel first, then real/integer columns only
    zDensCombo->addItem("-- (none) --");
    for (const auto& c : plottableX)
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
    for (const auto& c : plottableX)
        histCombo->addItem(colLabel(c));

    auto* histNBins = new QSpinBox;
    histNBins->setRange(1, 10000);
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

    // --- clearCanvas lambda ---
    // Removes the QCPColorScale (column 1 of plotLayout), all plottables, and all
    // annotation items (contour labels).  Resets the color cycling counters.
    auto clearCanvas = [=]() {
        if (customPlot->plotLayout()->hasElement(0, 1)) {
            customPlot->plotLayout()->remove(customPlot->plotLayout()->element(0, 1));
            customPlot->plotLayout()->simplify();
        }
        customPlot->clearPlottables();
        customPlot->clearItems();
        *hasContent      = false;
        *contourColorIdx = 0;
        *scatterColorIdx = 0;
        customPlot->replot();
    };

    // --- doScatter lambda ---
    auto doScatter = [=]() {
        const int xIdx = xScatterCombo->currentIndex();
        const bool useRowIdx = (xIdx == 0);
        const std::string xKey = useRowIdx
            ? std::string()
            : plottableX[static_cast<size_t>(xIdx - 1)].name;
        const std::string yKey =
            plottableY[static_cast<size_t>(yScatterCombo->currentIndex())].name;

        std::vector<std::shared_ptr<QCPGraph>> graphs;
        try {
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
            g->setPen(QPen(col));
        }

        const QString xLabel = useRowIdx
            ? QString("Row index")
            : colLabel(plottableX[static_cast<size_t>(xIdx - 1)]);
        const QString yLabel = colLabel(plottableY[static_cast<size_t>(yScatterCombo->currentIndex())]);
        customPlot->xAxis->setLabel(xLabel);
        customPlot->yAxis->setLabel(yLabel);

        customPlot->rescaleAxes();
        customPlot->replot(QCustomPlot::rpQueuedReplot);
        *hasContent = true;
    };

    // --- doHistogram lambda ---
    auto doHistogram = [=]() {
        if (plottableX.empty()) {
            QMessageBox::warning(nullptr, "Histogram", "No suitable columns available.");
            return;
        }
        const std::string colName =
            plottableX[static_cast<size_t>(histCombo->currentIndex())].name;

        const size_t nBins   = static_cast<size_t>(histNBins->value()); // 0 = auto
        const bool   logBins = histLogX->isChecked();
        const bool   logY    = histLogY->isChecked();

        std::vector<std::shared_ptr<QCPGraph>> graphs;
        try {
            graphs = QtFITSviewer::QtHistogram(customPlot, tab, colName, nBins);
        } catch (const std::exception& e) {
            QMessageBox::critical(nullptr, "Histogram error", QString::fromStdString(e.what()));
            return;
        }

        // Apply log scale types after the fact (QtHistogram's FITStable overload
        // hardcodes logFlag=0; we set the axis scale directly on the plot).
        customPlot->xAxis->setScaleType(logBins ? QCPAxis::stLogarithmic : QCPAxis::stLinear);
        customPlot->yAxis->setScaleType(logY    ? QCPAxis::stLogarithmic : QCPAxis::stLinear);

        const QColor col = scatterPalette[(*scatterColorIdx)++ % nScatterColors];
        for (auto& g : graphs)
            g->setPen(QPen(col));

        customPlot->xAxis->setLabel(colLabel(plottableX[static_cast<size_t>(histCombo->currentIndex())]));
        customPlot->yAxis->setLabel("count");

        customPlot->rescaleAxes();
        customPlot->replot(QCustomPlot::rpQueuedReplot);
        *hasContent = true;
    };

    // --- doDensity lambda ---
    // When Retain is unchecked: draws a filled QCPColorMap (QtDensityPlot).
    // When Retain is checked: computes the density grid, then discards the color map
    // and overlays QCustomContours contour lines instead, using a new uniform color
    // on each call so successive overlaid contour layers are visually distinct.
    auto doDensity = [=](bool retain) {
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
                ? plottableY[static_cast<size_t>(yDensIdx - 1)].name
                : plottableY[static_cast<size_t>(xDensIdx - 1)].name;
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
                    ? plottableY[static_cast<size_t>(yDensIdx - 1)].name
                    : plottableY[static_cast<size_t>(xDensIdx - 1)].name;
                colorMap = QtFITSviewer::QtDensityPlot(customPlot, tab, singleCol,
                                                        nBinsX, nBinsY, logFlag);
            } else {
                const std::string xCol = plottableY[static_cast<size_t>(xDensIdx - 1)].name;
                const std::string yCol = plottableY[static_cast<size_t>(yDensIdx - 1)].name;
                const int zIdx = zDensCombo->currentIndex();
                const std::string zCol = (zIdx == 0)
                    ? std::string()
                    : plottableX[static_cast<size_t>(zIdx - 1)].name;
                colorMap = QtFITSviewer::QtDensityPlot(customPlot, tab, xCol, yCol, zCol,
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

        if (!retain) {
            // No-retain path: keep the filled color map as-is.
            // Axis labels and replot were already handled inside QtDensityPlot.
            *hasContent = true;
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
        *hasContent = true;
    };

    // --- Right panel: tab widget + Plot/Clear/Retain below it ---
    auto* rightPanel  = new QWidget;
    auto* rightLayout = new QVBoxLayout(rightPanel);
    rightLayout->setContentsMargins(4, 4, 4, 4);
    rightLayout->addWidget(tabWidget);

    auto* plotBtn  = new QPushButton("Plot");
    auto* clearBtn = new QPushButton("Clear");
    auto* retainBox = new QCheckBox("Retain plot");

    auto* btnLayout = new QHBoxLayout;
    btnLayout->addWidget(plotBtn);
    btnLayout->addWidget(clearBtn);
    rightLayout->addLayout(btnLayout);
    rightLayout->addWidget(retainBox);
    rightPanel->setMinimumWidth(160);

    // --- Central container: splitter only ---
    auto* splitter = new QSplitter(Qt::Horizontal);
    splitter->addWidget(leftPanel);
    splitter->addWidget(customPlot);
    splitter->addWidget(rightPanel);
    splitter->setSizes({100, 600, 170});

    win->setCentralWidget(splitter);
    win->resize(1100, 680);

    // --- Wire connections ---
    // Plot button: dispatch to the active tab's draw function
    QObject::connect(plotBtn, &QPushButton::clicked, win, [=]() {
        const bool retain = retainBox->isChecked();
        if (!retain) clearCanvas();
        switch (tabWidget->currentIndex()) {
            case 0: doScatter();      break;
            case 1: doDensity(retain); break;
            case 2: doHistogram();    break;
            default: break;
        }
    });

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
        FITSmanager fm(path);
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
    FITSmanager fm(file);
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
                FITSmanager fm(m_currentFile);
                std::shared_ptr<FITScube> cube = fm.GetImageAtIndex(hduIndex);
                if (cube->GetDimension() < 3)
                    throw std::runtime_error("Selected HDU is not a 3-D cube.");
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

            FITSmanager ffile(in_file[0]);
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

        FITSmanager ffile(in_file[0]);
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