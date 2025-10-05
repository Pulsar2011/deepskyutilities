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

#endif

#include <DSTCore/OptionManager.h>
#include <DSTfits/FITSmanager.h>
#include <DSTfits/FITSimg.h>
#include <DSTfits/FITShdu.h>


#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future> // added for std::async

typedef std::thread THREAD;

std::thread _draw_histo;
std::thread _print_;

using namespace DSL;

typedef const std::vector<std::string> flist;
typedef std::map<char,std::pair<size_t,size_t> > roi;

size_t histo_bins=100;

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
        AddOption("--hdu", "int","Header block to be open. If not given, the primary header will be opened");
        AddOption("--opt/++opt","sqrtz/logz/quadz","Define the color ramp used to display the image.");
        AddOption("--bins/-b","int","Number of bin to display the histograms of the images.");
        AddOption("--roi","x0 xsize y0 ysize","Define a refgion of interest to be ploted. Other data outside this ROI are disregarded. The ROI shall be given in pixel not in WCS");
        AddOption("--print","Print pixel content to stdout.");
        AddOption("--info","Print statistical information on the array distribution.");
        AddOption("--draw","Disable drawing!");
        AddParam("Fits_Filename","FITS file URL and name that will be opened for analyses.");    
        
        DOption::prg_name += "FitsViewer";
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


void mouseHasMove(QCustomPlot *plt,const std::shared_ptr<const FITScube>& in, QCPGraph *x,QCPGraph *y)
{
    
    std::vector<double> px_coo({plt->xAxis->pixelToCoord(plt->mapFromGlobal(QCursor::pos()).x()),plt->yAxis->pixelToCoord(plt->mapFromGlobal(QCursor::pos()).y())});

    if (px_coo[0] < 0 || px_coo[0] >= in->Size(1) ||
        px_coo[1] < 0 || px_coo[1] >= in->Size(2) )
        return;

    try
    {
        QtFITSviewer::QtProfileX(plt, in,static_cast<size_t>(px_coo[1]+0.5),x);
        QtFITSviewer::QtProfileY(plt, in,static_cast<size_t>(px_coo[0]+0.5),y, true);
    }
    catch(std::out_of_range& e)
    {
        if((DSL::verbose&DSL::verboseLevel::VERBOSE_DEBUG) == DSL::verboseLevel::VERBOSE_DEBUG)
        std::cerr << "[WARNING] mouseHasMove: " << e.what() << std::endl;
    }

}

QCustomPlot* Draw(const std::shared_ptr<const FITScube>& in, std::string opt)
{
    if(UOption::ddraw == UOption::drawing_mode::OFF)
        return NULL;
    
    QtFITSviewer *mw = new QtFITSviewer(400, 250, 1200, 900);
    mw->setWindowTitle(in->GetName().c_str());

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
    
    image->setMaximumSize(850, 850);
    image->setMinimumSize(100, 100);
    yproj->setMinimumSize(10 , 100);
    xproj->setMinimumSize(100, 10);
    yproj->setMaximumSize(100, 850);
    xproj->setMaximumSize(850, 100);
    histo->setMaximumSize(850, 850);
    histo->setMinimumSize(50, 50);
    
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
    

    _draw_histo= std::thread([=]() {QtFITSviewer::QtHistogram(customPlot, in, static_cast<size_t>(histo_bins), px_pdf, px_cdf);});
    

    customPlot->connect(image->axis(QCPAxis::atBottom), SIGNAL(rangeChanged(QCPRange)), xproj->axis(QCPAxis::atBottom), SLOT(setRange(QCPRange)));
    
    customPlot->connect(image->axis(QCPAxis::atLeft), SIGNAL(rangeChanged(QCPRange)), yproj->axis(QCPAxis::atRight), SLOT(setRange(QCPRange)));
    
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes | QCP::iSelectOther);
    //customPlot->replot();
    
    auto f = std::bind(mouseHasMove,customPlot,in,x_projection,y_projection);
    QAbstractAnimation::connect(customPlot, &QCustomPlot::mouseMove, f);
    
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
    
        if(in_file.size() < 2)
            throw std::runtime_error("Not enougth number of argument. Please provide filename to open.");
    
        if(in_file.size() > 2)
            throw std::runtime_error("Too many argument! I can handle more than 1 file.");
    
        //3- OPEN IMAGE AND HEADER BLOCK
        std::shared_ptr<FITScube> img   = nullptr;
    
        FITSmanager ffile(in_file[1]);
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
        
        if(plt!=NULL)
        {        
            {
                std::lock_guard<std::mutex> lk(QtFITS::io_mutex);
                foreach (QWidget *widget, QApplication::allWidgets())
                {
                    widget->setGeometry(widget->x(), widget->y(), widget->width()*1.25, widget->height()/((img->Size(2)>=500)?1.25:1));
                    widget->update();
                    plt->replot();
                }
            }
        }

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
            _draw_histo.join();
        
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


