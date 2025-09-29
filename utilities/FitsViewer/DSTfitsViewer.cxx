#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <map>

#ifdef __WITH_QT__

#include <QtConcurrent/QtConcurrent>
#include <DSTdisplay/QtFitsPloting.h>
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

std::mutex io_mutex;
std::mutex msg_m;
std::condition_variable msg_cv;
int finished = -1;

using namespace DSL;

typedef const std::vector<std::string> flist;
typedef std::map<char,std::pair<size_t,size_t> > roi;

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
    size_t NBINS;
    std::string opt_string;
    rmode MODE;
    bool do_Print;
    double scale;
    bool doScale;
    bool a_info;

public:
    
    static drawing_mode ddraw;
    roi ROI;
    
    UOption():DOption(),HDU(1),RBIN(0),NBINS(100),opt_string(),MODE(rmode::NONE), do_Print(false),scale(1.),doScale(false),a_info(false),ROI()
    {
        AddDescription("Open a FITS header block and display data.");
        AddOption("--hdu", "int","Header block to be open. If not given, the primary header will be opened");
        AddOption("--opt/++opt","sqrtz/logz/quadz","Define the color ramp used to display the image.");
        AddOption("--bins/-b","int","Number of bin to display the histograms of the images.");
        AddOption("--roi","x0 xsize y0 ysize","Define a refgion of interest to be ploted. Other data outside this ROI are disregarded. The ROI shall be given in pixel not in WCS");
        AddOption("--print","Print pixel content to stdout.");
        AddOption("--info","Print statistical information on the array distribution.");
        AddOption("--draw","Disable drawing!");
        AddOption("--scale/-s","scaling factor","Normalisation factor to be used on the images.");
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
            NBINS = std::stoi(*it);
        }
        else if(!it->compare("--scale") || !it->compare("++scale") || !it->compare("-s") || !it->compare("+s"))
        {
            if(std::next(it) != arg.cend()){it++;}else{return false;}
            scale = std::stod(*it);
            doScale |= true;
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
    inline const size_t&       N_bins() const {return NBINS;}
    inline const rmode&        Mode()  const {return MODE;}
    inline const int16_t&      Header()const {return HDU;}
    inline const std::string&  draw()  const {return opt_string;}
    inline bool doPrint() const {return do_Print;}
    inline bool ArrayInfo() const {return a_info;}
    inline const double& Scale() const {return scale;}
    inline const bool& DoScale() const {return doScale;}
    
};

UOption::drawing_mode UOption::ddraw = UOption::drawing_mode::LINEAR;

#ifdef __WITH_QT__

template<typename T>
QSharedPointer< QCPGraphDataContainer > get_x_profile(const std::shared_ptr<const FITScube>& in, const std::slice& s, const std::vector<size_t>& px,std::vector<double>& min_max)
{
    min_max.clear(); while(min_max.size() != 2) min_max.push_back(0);
    
    QCPGraphDataContainer* data = new QCPGraphDataContainer();
    const std::valarray<T>* _val_ = in->GetData<T>();
    std::valarray<T> val = (*_val_)[s];
    
    switch (UOption::ddraw)
    {
        case UOption::drawing_mode::LOG:
            val[val < 0.] = 0;
            val = log10(val);
            break;
            
        case UOption::drawing_mode::QUAD:
            val *= val;
            break;
            
        case UOption::drawing_mode::SQRT:
            val[val < 0.] = 0;
            val = sqrt(val);
            break;
            
        default:
            break;
    }
    
    for(size_t ix = 0; ix < val.size(); ix++)
    {
        double X= in->WorldCoordinates(std::vector<size_t>({ix+px[0],px[1]}))[0];
        double Y= val[ix];

        min_max[0] = (Y < min_max[0])?Y:min_max[0];
        min_max[1] = (Y > min_max[1])?Y:min_max[1];
        data->add(QCPGraphData(X,Y));
    }
    
    if(UOption::ddraw == UOption::drawing_mode::LOG)
        min_max[0] = static_cast<double>(((std::valarray<T>)val[val>0.]).min())*(1-0.5);
    
    return QSharedPointer< QCPGraphDataContainer >(data);
}

template<typename T>
QSharedPointer< QCPGraphDataContainer > get_y_profile(const std::shared_ptr<const FITScube>& in, const std::slice& s, const std::vector<size_t>& px,std::vector<double>& min_max)
{
    min_max.clear(); while(min_max.size() != 2) min_max.push_back(0);
    
    if(s.start()+s.stride()*s.size() >= in->Size())
        throw std::out_of_range("get_y_profile : stride is out of range");
    
    QCPGraphDataContainer* data = new QCPGraphDataContainer();
    const std::valarray<T>* _val_ = in->GetData<T>();
    std::valarray<T> val = (*_val_)[s];
    
    switch (UOption::ddraw)
    {
        case UOption::drawing_mode::LOG:
            val[val < 0.] = 0;
            val = log10(val);
            break;
            
        case UOption::drawing_mode::QUAD:
            val *= val;
            break;
            
        case UOption::drawing_mode::SQRT:
            val[val < 0.] = 0;
            val = sqrt(val);
            break;
            
        default:
            break;
    }
    
    for(size_t ix = 0; ix < val.size(); ix++)
    {
        double X= in->WorldCoordinates(std::vector<size_t>({px[0],ix+px[1]}))[1];
        double Y= val[ix];
        
        min_max[0] = (Y < min_max[0])?Y:min_max[0];
        min_max[1] = (Y > min_max[1])?Y:min_max[1];
        data->add(QCPGraphData(X,Y));
    }
    
    if(UOption::ddraw == UOption::drawing_mode::LOG)
        min_max[0] = static_cast<double>(((std::valarray<T>)val[val>0.]).min())*(1-0.5);
    
    return QSharedPointer< QCPGraphDataContainer >(data);
}

void mouseHasMove(QCustomPlot *plt,const std::shared_ptr<const FITScube>& in, QCPGraph *x,QCPGraph *y)
{
    std::vector<double> WCS(2);
    std::vector<double> px_coo;
    
    size_t x_start, x_stop;
    size_t y_start, y_stop;
  
    size_t dx, dy;
    
    WCS[0] = plt->xAxis->pixelToCoord(plt->mapFromGlobal(QCursor::pos()).x());
    WCS[1] = plt->yAxis->pixelToCoord(plt->mapFromGlobal(QCursor::pos()).y());
    
    if(WCS[0] < plt->xAxis->range().lower || WCS[0] > plt->xAxis->range().upper)
    return;
    
    if(WCS[1] < plt->yAxis->range().lower || WCS[1] > plt->yAxis->range().upper)
    return;
    
    px_coo  = in->World2Pixel(WCS);
    size_t px = in->PixelIndex(px_coo);

    plt->window()->setStatusTip(QString("[")+QString::number(static_cast<size_t>(px_coo[0]))+QString(";")+QString::number(static_cast<size_t>(px_coo[1]))+QString("] : (")+QString::number(WCS[0])+QString(";")+QString::number(WCS[1])+QString(") I: ")+QString::number(in->FloatValueAtPixel(px)));
    //plt->statusBar()->showMessage("");

    double img_line_start = in->WorldCoordinates(std::vector<double>({0,px_coo[1]}))[0];
    double img_line_stop  = in->WorldCoordinates(std::vector<double>({static_cast<double>(in->Size(1)-1),px_coo[1]}))[0];
    if(img_line_stop < img_line_start)
        std::swap(img_line_start,img_line_stop);
    
    double img_col_start = in->WorldCoordinates(std::vector<double>({px_coo[0],0}))[1];
    double img_col_stop  = in->WorldCoordinates(std::vector<double>({px_coo[0],static_cast<double>(in->Size(2)-1)}))[1];
    if(img_col_stop < img_col_start)
        std::swap(img_col_start,img_col_stop);
    
    x_start = in->PixelIndex(in->World2Pixel({std::max(plt->xAxis->range().lower,img_line_start),WCS[1]}));
    y_start = in->PixelIndex(in->World2Pixel({WCS[0],std::max(plt->yAxis->range().lower,img_col_start)}));
    
    x_stop = in->PixelIndex(in->World2Pixel({std::min(plt->xAxis->range().upper       ,img_line_stop),WCS[1]}));
    y_stop = in->PixelIndex(in->World2Pixel({WCS[0],img_col_stop}));
    
    if(x_start > x_stop)
        std::swap(x_start,x_stop);
    
    if(y_start > y_stop)
        std::swap(y_start,y_stop);
    
    dx      = static_cast<size_t>(std::min(std::abs(static_cast<double>(x_stop)-static_cast<double>(x_start)),static_cast<double>(in->Size(1)-1)));
    dy      = static_cast<size_t>(std::min(std::abs(static_cast<double>(y_stop)-static_cast<double>(y_start)),static_cast<double>(in->Size(2)-1)));
    
    std::vector<size_t> coo_x = in->PixelCoordinates(x_start);
    std::vector<size_t> coo_y = in->PixelCoordinates(y_start);
    
    dx = static_cast<size_t>((static_cast<double>(in->Size(1))-1)-coo_x[0]);
    dy = static_cast<size_t>((static_cast<double>(in->Size(2))-1)-coo_y[1]);
    
    std::vector<double> x_range;
    std::vector<double> y_range;
    try
    {
        int eqBitPix= in->GetEqBitPerPixel();
        switch (eqBitPix)
        {
            case BYTE_IMG:
                x->setData(get_x_profile<uint8_t>(in,std::slice(x_start,dx,1)          ,coo_x, x_range));
                y->setData(get_y_profile<uint8_t>(in,std::slice(y_start,dy,in->Size(1)),coo_y, y_range));
                break;

            case SBYTE_IMG:
                x->setData(get_x_profile<int8_t>(in,std::slice(x_start,dx,1)          ,coo_x, x_range));
                y->setData(get_y_profile<int8_t>(in,std::slice(y_start,dy,in->Size(1)),coo_y, y_range));
                break;
            
            case USHORT_IMG:
                x->setData(get_x_profile<uint16_t>(in,std::slice(x_start,dx,1)          ,coo_x, x_range));
                y->setData(get_y_profile<uint16_t>(in,std::slice(y_start,dy,in->Size(1)),coo_y, y_range));
                break;

            case SHORT_IMG:
                x->setData(get_x_profile<int16_t>(in,std::slice(x_start,dx,1)          ,coo_x, x_range));
                y->setData(get_y_profile<int16_t>(in,std::slice(y_start,dy,in->Size(1)),coo_y, y_range));
                break;

            case ULONG_IMG:
                x->setData(get_x_profile<uint32_t>(in,std::slice(x_start,dx,1)          ,coo_x, x_range));
                y->setData(get_y_profile<uint32_t>(in,std::slice(y_start,dy,in->Size(1)),coo_y, y_range));
                break;

            case LONG_IMG:
                x->setData(get_x_profile<int32_t>(in,std::slice(x_start,dx,1)          ,coo_x, x_range));
                y->setData(get_y_profile<int32_t>(in,std::slice(y_start,dy,in->Size(1)),coo_y, y_range));
                break;

            case ULONGLONG_IMG:
                x->setData(get_x_profile<uint64_t>(in,std::slice(x_start,dx,1)          ,coo_x, x_range));
                y->setData(get_y_profile<uint64_t>(in,std::slice(y_start,dy,in->Size(1)),coo_y, y_range));
                break;

            case LONGLONG_IMG:
                x->setData(get_x_profile<int64_t>(in,std::slice(x_start,dx,1)          ,coo_x, x_range));
                y->setData(get_y_profile<int64_t>(in,std::slice(y_start,dy,in->Size(1)),coo_y, y_range));
                break;

            case FLOAT_IMG:
                x->setData(get_x_profile<float>(in,std::slice(x_start,dx,1)          ,coo_x, x_range));
                y->setData(get_y_profile<float>(in,std::slice(y_start,dy,in->Size(1)),coo_y, y_range));
                break;

            case DOUBLE_IMG:
                x->setData(get_x_profile<double>(in,std::slice(x_start,dx,1)          ,coo_x, x_range));
                y->setData(get_y_profile<double>(in,std::slice(y_start,dy,in->Size(1)),coo_y, y_range));
                break;
        
            default:
                break;
        }

    }
    catch(...)
    {
        return;
    }
    
    x->valueAxis()->setRange(x_range[0],x_range[1]);
    y->valueAxis()->setRange(y_range[1],y_range[0]);
    
    plt->replot();
    
    WCS.clear();
    coo_x.clear();
    coo_y.clear();
}

template<typename T>
void DrawHistogram(const std::shared_ptr<const FITScube>& in, size_t n_bins, QCPGraph *px_pdf, QCPGraph *px_cdf, QCustomPlot *customPlot)
{
    double max = in->GetMaximum();
    double min = in->GetMinimum();
    
    switch (UOption::ddraw)
    {
        case UOption::drawing_mode::LOG:
            max=log10(max);
            min=(min > 0)?log10(min):0;
            break;
            
        case UOption::drawing_mode::QUAD:
            max *= max;
            min *= min;
            break;
            
        case UOption::drawing_mode::SQRT:
            max=sqrt(max);
            min=(min > 0)?sqrt(min):0;
            break;
            
        default:
            break;
    }
    
    double SPAN = min/max;
    
    QVector<double> val, entry, sum;
    
    if(SPAN > 1000)
    {
        for(size_t i = 0; i <= n_bins; i++)
            val.push_back( pow(10.,static_cast<double>(i)/(static_cast<double>(n_bins)-1)*log10(max/min)+log10(min)));
    }
    else
    {
        for(size_t i = 0; i <= n_bins; i++)
            val.push_back( static_cast<double>(i)/(static_cast<double>(n_bins)-1)*(max-min)+min);
    }
    
    const std::valarray<T>* _arr_ = in->GetData<T>();
    std::vector<T> array;
    array.assign(std::begin(*_arr_),std::end(*_arr_));
    std::sort(array.begin(),array.end());
    
    //switch (UOption::ddraw)
    //{
    //    case UOption::drawing_mode::LOG:
    //        array[array < 0.] = 0;
    //        array = std::log10(array);
    //        break;
    //        
    //    case UOption::drawing_mode::QUAD:
    //        array *= array;
    //        break;
    //        
    //    case UOption::drawing_mode::SQRT:
    //        array[array < 0.] = 0;
    //        array = std::sqrt(array);
    //        break;
    //        
    //    default:
    //        break;
    //}
    
    double total = static_cast<double>(array.size());
    double iMax  = 0;
    double this_sum   = 0;
    
    for(QVector<double>::const_iterator iBin = val.cbegin(); iBin != val.cend(); iBin++)
    {
        typename std::vector<T>::const_iterator i = array.cbegin();
        double this_val  = (*iBin);
        double this_entry = 0.;

        while(this_val < (*(iBin+1)) && i != array.cend())
        {
            switch (UOption::ddraw)
            {
                case UOption::drawing_mode::LOG:
                    this_val=( (*i) > 0 )?log10(*i):0;
                    break;
                case UOption::drawing_mode::QUAD:
                    this_val= (*i)*(*i);
                    break;
                case UOption::drawing_mode::SQRT:
                    this_val=( (*i) > 0 )?sqrt(*i):0;
                    break;
                default:
                    this_val=(*i);
                    break;
            }

            if(this_val < (*(iBin+1)))
            {
                this_sum++;
                this_entry++;
            }

            i++;
        }

        entry.push_back(this_entry);
        sum.push_back(this_sum);
        iMax = (this_entry>iMax)?this_entry:iMax;
    }
    
    //sum  .push_back(static_cast<double>(((std::valarray<T>)array[array <= static_cast<T>(max)]).size()));
    //entry.push_back(static_cast<double>(((std::valarray<T>)array[array >= static_cast<T>(max)]).size()));
    
    auto xVals   = val;
    auto pdfVals = entry;
    auto cdfVals = sum;

    QMetaObject::invokeMethod(customPlot, [px_pdf, px_cdf, customPlot,
                                           xVals   = std::move(xVals),
                                           pdfVals = std::move(pdfVals),
                                           cdfVals = std::move(cdfVals)]() mutable {
        std::lock_guard<std::mutex> lk(io_mutex);
        px_pdf->setData(xVals, pdfVals);
        px_cdf->setData(xVals, cdfVals);
        px_pdf->valueAxis()->rescale();
        px_pdf->keyAxis()->rescale();
        customPlot->replot();
        {
            std::lock_guard<std::mutex> dk(msg_m);
            if(finished == -1)
                finished = 1;
        }
        msg_cv.notify_one();
    }, Qt::QueuedConnection);
}

QCustomPlot* Draw(const std::shared_ptr<const FITScube>& in, std::string opt)
{
    if(UOption::ddraw == UOption::drawing_mode::OFF)
        return NULL;
    
    QtFITSviewer *mw = new QtFITSviewer();
    mw->setWindowTitle(in->GetName().c_str());

    //size_t xsize = std::min(static_cast<size_t>(900),in->Size(1)+200);
    //size_t ysize = std::min(static_cast<size_t>(900),in->Size(2)+200);

    mw->setGeometry(400, 250, 1200, 900);

    QCustomPlot *customPlot = mw->GetPlot();
    
    QtDisplay(customPlot,in,opt);
    
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
    
    int eqBits = in->GetEqBitPerPixel();
    switch (eqBits)
    {
        case BYTE_IMG:
        _draw_histo= std::thread(DrawHistogram<uint8_t>, in,100, px_pdf, px_cdf, customPlot);
        break;

        case SBYTE_IMG:
        _draw_histo= std::thread(DrawHistogram<int8_t>, in,100, px_pdf, px_cdf, customPlot);
        break;

        case SHORT_IMG:
        _draw_histo= std::thread(DrawHistogram<int16_t>, in,100, px_pdf, px_cdf, customPlot);
        break;

        case USHORT_IMG:
        _draw_histo= std::thread(DrawHistogram<uint16_t>, in,100, px_pdf, px_cdf, customPlot);
        break;

        case LONG_IMG:
        _draw_histo= std::thread(DrawHistogram<int32_t>, in,100, px_pdf, px_cdf, customPlot);
        break;

        case ULONG_IMG:
        _draw_histo= std::thread(DrawHistogram<uint32_t>, in,100, px_pdf, px_cdf, customPlot);
        break;

        case LONGLONG_IMG:
        _draw_histo= std::thread(DrawHistogram<int64_t>, in,100, px_pdf, px_cdf, customPlot);
        break;

        case ULONGLONG_IMG:
        _draw_histo= std::thread(DrawHistogram<uint64_t>, in,100, px_pdf, px_cdf, customPlot);
        break;

        case FLOAT_IMG:
        _draw_histo= std::thread(DrawHistogram<float>, in,100, px_pdf, px_cdf, customPlot);
        break;

        case DOUBLE_IMG:
        _draw_histo= std::thread(DrawHistogram<double>, in,100, px_pdf, px_cdf, customPlot);
        break;
    
        default:
        break;
    }

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
        
        if(fOPT.DoScale())
        {
            (*img)/=fOPT.Scale();
        }
        
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
            msg_m.lock();
            std::cout<<"\033[01m=====\033[34m"<< fOPT.GetFiles()[1]<<"\033[0m"<<std::endl
                   <<"\033[31m     |- AVG :\033[0m"<<mean<<" "<<((img->HDU().Exists("BUNIT"))?img->HDU().GetValueForKey("BUNIT"):"")<<std::endl
                   <<"\033[31m     |- MED :\033[0m"<<median<<" "<<((img->HDU().Exists("BUNIT"))?img->HDU().GetValueForKey("BUNIT"):"")<<std::endl
                   <<"\033[31m     |- 25% :\033[0m"<<pp25<<" "<<((img->HDU().Exists("BUNIT"))?img->HDU().GetValueForKey("BUNIT"):"")<<std::endl
                   <<"\033[31m     |- 75% :\033[0m"<<pp75<<" "<<((img->HDU().Exists("BUNIT"))?img->HDU().GetValueForKey("BUNIT"):"")<<std::endl;
            msg_m.unlock();
        }

        
        if( fOPT.ArrayInfo() )
        {
            msg_m.lock();
            std::cout<<"\033[32m==== "<<in_file[0]<<"["<<fOPT.Header()<<"]\033[0m"<<std::endl;
            std::cout<<" FILE |  NAME | MEAN | MEDIAN | Q1 | Q3 | RMS | MIN | MAX | SUM"<<std::endl;
            std::cout<<fOPT.GetFiles()[0]<<" | "<<img->GetName()<<" | "<<mean<<" | "<<median<<" | "<<pp25<<" | "<<pp75<<" | "<<rms<<" | "<<min<<" | "<<max<<" | "<<sum<<std::endl;
            msg_m.unlock();
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
                std::unique_lock<std::mutex> lk(msg_m);
                msg_cv.wait(lk, []{return finished != -1;});
            }
        
            {
                std::lock_guard<std::mutex> lk(io_mutex);
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


