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
#include <QFileDialog>
#include <DSTdisplay/QCustomContours.h>

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
        if(DSD::gGIT::version() < "v1.1.0")
            throw std::runtime_error("[ERROR] UOption::UOption: Incompatible DSD version, please update your DSD installation to v1.1.0 at minimum.");
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

        if (in_file.size() > 1)
            throw std::runtime_error("Too many arguments! I can handle only 1 file.");

        // Headless path: --draw disables GUI
        if (UOption::ddraw == UOption::drawing_mode::OFF)
        {
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