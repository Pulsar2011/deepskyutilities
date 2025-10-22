#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <map>

#ifdef __WITH_QT__

#include <QtConcurrent/QtConcurrent>
#include <DSTdisplay/QtFITSviewer.h>
#include <QMetaObject>
#include <QMainWindow> 
#include <QStatusBar>  

#endif

#include <DSTCore/OptionManager.h>
#include <DSTCore/DSS_version.h>
#include <DSTools/DST_version.h>
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
std::thread _print_info;
std::thread _print_;

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

public:
    
    static drawing_mode ddraw;
    roi ROI;
    
    UOption():DOption(),HDU(1),RBIN(0),opt_string(),MODE(rmode::NONE), do_Print(false),a_info(false),ROI()
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
            throw std::runtime_error("[ERROR] UOption::UOption: Incompatible DSS version, please update your DSS installation.");
    }
    
    virtual ~UOption(){};
    
    bool UserParameters(std::vector<std::string>::const_iterator& it, const std::vector<std::string>& arg)
    {
        
        if(!it->compare("--hdu"))
        {
            if(std::next(it) != arg.cend()){it++;}else{return false;}
            HDU = std::stoi(*it);
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
    
    std::vector<double> px_coo({plt->xAxis->pixelToCoord(plt->mapFromGlobal(QCursor::pos()).x()),plt->yAxis->pixelToCoord(plt->mapFromGlobal(QCursor::pos()).y())});

    if (px_coo[0] < 0 || px_coo[0] >= in->Size(1) ||
        px_coo[1] < 0 || px_coo[1] >= in->Size(2) )
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

        // find top-level window (the QtFITSviewer is a QMainWindow) and write to its status bar
        QWidget *top = plt->window();
        if(top)
        {
            QMainWindow *mainWin = qobject_cast<QMainWindow*>(top);
            if(mainWin && mainWin->statusBar())
            {
                QString msg = QString("x=%1  y=%2  val=%3").arg(ix).arg(iy).arg(valStr);
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
                                                    ell->setMaximumSize(500, 1200);
                                                    

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
    customPlot->plotLayout()->addElement(1, 0, xproj_layout); // insert axis rect in first row
    
    QCPLayoutGrid *histo_layout = new QCPLayoutGrid;
    customPlot->plotLayout()->addElement(0, 3, histo_layout); // insert axis rect in first row
    
    QCPLayoutGrid *yproj_layout = new QCPLayoutGrid;
    customPlot->plotLayout()->addElement(0, 2, yproj_layout); // insert axis rect in first row
    
    
    // prepare axis rects that will be placed in the sublayout:
    QCPAxisRect *yproj = new QCPAxisRect(customPlot, false);
    QCPAxisRect *xproj = new QCPAxisRect(customPlot, false);
    QCPAxisRect *histo = new QCPAxisRect(customPlot, false);
    
    image->setMinimumSize(200, 200);
    yproj->setMinimumSize(50, 200);
    yproj->setMaximumSize(100, 1200);
    xproj->setMinimumSize(200, 50);
    xproj->setMaximumSize(1200, 100);
    histo->setMaximumSize(500, 1200);
    
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

    image->axis(QCPAxis::atBottom)->setScaleRatio(image->axis(QCPAxis::atLeft),1.0);
    
    mw->SetPlot(customPlot);
    mw->show();
    return customPlot;
}

#endif

int main(int argc, char **argv)
{
#ifdef __WITH_QT__
    
    try
    {
        //-1 READ OPTIONAL PARAMETERS
        UOption fOPT = UOption();
        if(!fOPT.InitParameters(argc, argv))
            return EXIT_FAILURE;
        
        std::ostream* fp = &std::cout;
        std::ofstream fout;
        if (fOPT.GetOutput().size() > 1)
        {
            fout.open(fOPT.GetOutput());
            fp = &fout;
        }
    
        //2- GET IMAGE TO BE READS
        int16_t debugLevel = fOPT.DebugLevel();
        if(debugLevel < 1)
            DSL::verbose = DSL::verboseLevel::VERBOSE_NONE;
        else if(debugLevel == 1)
            DSL::verbose = DSL::verboseLevel::VERBOSE_BASIC;
        else if(debugLevel == 2)
            DSL::verbose = DSL::verboseLevel::VERBOSE_DETAIL;
        else
            DSL::verbose = DSL::verboseLevel::VERBOSE_DEBUG;

        flist in_file = fOPT.GetFiles();
    
        if(in_file.size() < 1)
            throw std::runtime_error("Not enougth number of argument. Please provide filename to open.");
    
        if(in_file.size() > 1)
            throw std::runtime_error("Too many argument! I can handle more than 1 file.");
    
        //3- OPEN IMAGE AND HEADER BLOCK
        std::shared_ptr<FITScube> img   = nullptr;
    
        FITSmanager ffile(in_file[0]);
        img = ffile.GetImageAtIndex(fOPT.Header());
        
        if( fOPT.ROI.size() > 1)
        {
            size_t X0 = 0;
            size_t nX = img->Size(1);
            
            if(fOPT.ROI.find('x') != fOPT.ROI.end())
            {
                X0 = fOPT.ROI.find('x')->second.first;
                nX = fOPT.ROI.find('x')->second.second;
            }
            
            size_t Y0 = 0;
            size_t nY = img->Size(2);
            
            if(fOPT.ROI.find('y') != fOPT.ROI.end())
            {
                Y0 = fOPT.ROI.find('y')->second.first;
                nY = fOPT.ROI.find('y')->second.second;
            }
            
            std::shared_ptr<FITScube>   tmp = img->Window(X0, Y0, nX, nY);
            std::swap(img,tmp);
            tmp.reset();
        }

        // prepare variables and compute them in parallel
        double sum=0,mean=0, median=0, rms=0, min=0, max=0, pp25=0, pp75=0;

        auto fut_sum   = std::async(std::launch::async, [img]{ return img->GetSum(); });
        auto fut_mean  = std::async(std::launch::async, [img]{ return img->GetMean(); });
        auto fut_median= std::async(std::launch::async, [img]{ return img->GetMedian(); });
        auto fut_rms   = std::async(std::launch::async, [img]{ return img->GetVariance(); }); // consistent with prior "RMS" usage
        auto fut_min   = std::async(std::launch::async, [img]{ return img->GetMinimum(); });
        auto fut_max   = std::async(std::launch::async, [img]{ return img->GetMaximum(); });
        auto fut_pp25   = std::async(std::launch::async, [img]{ return img->Getpercentil(0.05); });
        auto fut_pp75  = std::async(std::launch::async, [img]{ return img->Getpercentil(0.95); });

        sum   = fut_sum.get();
        mean  = fut_mean.get();
        median= fut_median.get();
        rms   = fut_rms.get();
        min   = fut_min.get();
        max   = fut_max.get();
        pp25  = fut_pp25.get();
        pp75  = fut_pp75.get();

        if((DSL::verbose & DSL::verboseLevel::VERBOSE_BASIC)==DSL::verboseLevel::VERBOSE_BASIC)
        {
            QtFITS::msg_m.lock();
            std::cout<<"\033[01m=====\033[34m"<< fOPT.GetFiles()[1]<<"\033[0m"<<std::endl
                   <<"\033[31m     |- AVG :\033[0m"<<mean<<" "<<((img->HDU().Exists("BUNIT"))?img->HDU().GetValueForKey("BUNIT"):"")<<std::endl
                   <<"\033[31m     |- MED :\033[0m"<<median<<" "<<((img->HDU().Exists("BUNIT"))?img->HDU().GetValueForKey("BUNIT"):"")<<std::endl
                   <<"\033[31m     |- 25% :\033[0m"<<pp25<<" "<<((img->HDU().Exists("BUNIT"))?img->HDU().GetValueForKey("BUNIT"):"")<<std::endl
                   <<"\033[31m     |- 75% :\033[0m"<<pp75<<" "<<((img->HDU().Exists("BUNIT"))?img->HDU().GetValueForKey("BUNIT"):"")<<std::endl;
            QtFITS::msg_m.unlock();
        }

        
        if( fOPT.ArrayInfo() )
        {
            QtFITS::msg_m.lock();
            std::cout<<"\033[32m==== "<<in_file[0]<<"["<<fOPT.Header()<<"]\033[0m"<<std::endl;
            std::cout<<" FILE |  NAME | MEAN | MEDIAN | Q1 | Q3 | RMS | MIN | MAX | SUM"<<std::endl;
            std::cout<<fOPT.GetFiles()[0]<<" | "<<img->GetName()<<" | "<<mean<<" | "<<median<<" | "<<pp25<<" | "<<pp75<<" | "<<rms<<" | "<<min<<" | "<<max<<" | "<<sum<<std::endl;
            QtFITS::msg_m.unlock();
        }
        
        QApplication a(argc, argv);
        auto c_img = std::const_pointer_cast<const FITScube>(img);

        QCustomPlot* plt = Draw(c_img,fOPT.draw());
        
        foreach (QWidget *widget, QApplication::allWidgets())
        {
            if(plt == NULL)
                continue;
            widget->setGeometry(widget->x(), widget->y(), widget->width()*1.25, widget->height()/((img->Size(2)>=500)?1.25:1));
            widget->update();
            plt->replot();
        }
        
         if(plt!=NULL)
            a.exec();
        
//        if(plt!=NULL)
//        {        
//            {
//                std::lock_guard<std::mutex> lk(QtFITS::io_mutex);
//                foreach (QWidget *widget, QApplication::allWidgets())
//                {
//                    widget->setGeometry(widget->x(), widget->y(), widget->width()*1.25, widget->height()/((img->Size(2)>=500)?1.25:1));
//                    widget->update();
//                    plt->replot();
//                }
//            }
//        }

        bool printing_is_done=false;
        if(fOPT.doPrint())
        {
               switch (img->GetEqBitPerPixel())
                {
                    case BYTE_IMG:
                    {
                        auto F = std::bind(&DSL::FITSimg<uint8_t>::Print, dynamic_cast< DSL::FITSimg<uint8_t>* >(img.get()));
                        _print_=std::thread(F);
                        break;
                    }

                    case SBYTE_IMG:
                    {
                        auto F = std::bind(&DSL::FITSimg<int8_t>::Print, dynamic_cast< DSL::FITSimg<int8_t>* >(img.get()));
                        _print_=std::thread(F);
                        break;
                    }

                    case SHORT_IMG:
                    {
                        auto F = std::bind(&DSL::FITSimg<int16_t>::Print, dynamic_cast< DSL::FITSimg<int16_t>* >(img.get()));
                        _print_=std::thread(F);
                        break;
                    }

                    case USHORT_IMG:
                    {
                        auto F = std::bind(&DSL::FITSimg<uint16_t>::Print, dynamic_cast< DSL::FITSimg<uint16_t>* >(img.get()));
                        _print_=std::thread(F);
                        break;
                    }

                    case LONG_IMG:
                    {
                        auto F = std::bind(&DSL::FITSimg<int32_t>::Print, dynamic_cast< DSL::FITSimg<int32_t>* >(img.get()));
                        _print_=std::thread(F);
                        break;
                    }

                    case ULONG_IMG:
                    {
                        auto F = std::bind(&DSL::FITSimg<uint32_t>::Print, dynamic_cast< DSL::FITSimg<uint32_t>* >(img.get()));
                        _print_=std::thread(F);
                        break;
                    }

                    case LONGLONG_IMG:
                    {
                        auto F = std::bind(&DSL::FITSimg<int64_t>::Print, dynamic_cast< DSL::FITSimg<int64_t>* >(img.get()));
                        _print_=std::thread(F);
                        break;
                    }

                    case ULONGLONG_IMG:
                    {
                        auto F = std::bind(&DSL::FITSimg<uint64_t>::Print, dynamic_cast< DSL::FITSimg<uint64_t>* >(img.get()));
                        _print_=std::thread(F);
                        break;
                    }

                    case FLOAT_IMG:
                    {
                        auto F = std::bind(&DSL::FITSimg<float>::Print, dynamic_cast< DSL::FITSimg<float>* >(img.get()));
                        _print_=std::thread(F);
                        break;
                    }

                    case DOUBLE_IMG:
                    {
                        auto F = std::bind(&DSL::FITSimg<double>::Print, dynamic_cast< DSL::FITSimg<double>* >(img.get()));
                        _print_=std::thread(F);
                        break;
                    }

                    default:
                    {
                        if(fOPT.GetOutput().size() > 0)
                            std::cout<<" ... WRITTING TO "<<fOPT.GetOutput()<<" plz, wait!"<<std::endl;
                        img->Print();
                        printing_is_done=true;
                    }
                }
        }
        
        if(plt != NULL)
        {
            _draw_histo.join();
            _print_info.join();
        }
        
        if(fOPT.doPrint() && !printing_is_done)
        {
            if(fOPT.GetOutput().size() > 0)
                std::cout<<" ... WRITTING TO "<<fOPT.GetOutput()<<" plz, wait!"<<std::endl;
            _print_.join();
        }
        
        img.reset();
        ffile.Close();

    }
    catch(std::exception& e)
    {
        std::cerr<<e.what()<<std::endl;
        throw;
    }
    catch(...)
    {
        return EXIT_FAILURE;
    }
#else
    std::cerr<<"\033[43;34m*** [FitsViewer]: Errors *** "<<std::endl
        <<"THIS APPLICATION REQUIRE Qt5+ TO BE INSTALLED"<<std::endl
        <<"PLZ, INSTALL Qt5 or Qt6 AND RETRY\033[0m"<<std::endl;
    return EXIT_FAILURE;
#endif
    
    return EXIT_SUCCESS;
}