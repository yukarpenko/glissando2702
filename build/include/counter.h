/** \file counter.h
 * Part of GLISSANDO
 * 
 */


#ifndef _GL_COUNTER
  #define _GL_COUNTER  


/*******************************
 counting classes
*******************************/

//! simplest counting class
class counter {
private:
    Int_t count; //!< number of entries
    Double_t value; //!< sum of values counted
public:
//! constructor
    counter(){count=0;value=0;};
//! destructor
    ~counter(){};
//! reset the counter 
    void reset(){count=0;value=0;}; 
//! add entry to the counter
    void add(
            Double_t s //!< value added
            ){count++; value=value+s;};
//! get the number of entries
    Int_t getN(){return count;};
//! get the sum of values
    Double_t get(){return value;};
//! get the mean value
    Double_t mean(){return value/count;};
};

//! counting class with variance
class counter2 {
private:
  Int_t     count;  //!< number of entries
  Double_t  value,  //!< sum of values
            value2; //!< sum of squares of values
public:
//! constructor
    counter2(){count=0;value=0;value2=0;};
//! destructor
    ~counter2(){};
//! reset the counter 
    void reset(){count=0;value=0;value2=0;};
//! add entry
    void add(
            Double_t s //!< value added
            ){count++; value=value+s; value2=value2+s*s;};
//! get the number of entries
    Int_t getN(){return count;};
//! get the sum of values
    Double_t get(){return value;};
//! get the sum of squares of values
    Double_t get2(){return value2;};
//! get the mean value
    Double_t mean(){return value/count;};
//! get the variance
    Double_t var(){return value2/(count-1)-value/count*value/(count-1);}
//! get the variance multiplied with (N-1)/N
    Double_t vara(){return value2/count-value/count*value/count;};};

//! 2-dimensional counting class
class counter_2D {
private:
  Int_t count;          //!< number of entries
  Double_t     valuex,  //!< sum of x values
               valuey,  //!< sum of y values
               valuex2, //!< sum of x^2
               valuey2, //!< sum of y^2
               valuexy; //!< sum of x*y
public:
//! constructor
    counter_2D(){count=0;valuex=0;valuey=0;valuex2=0;valuey2=0;valuexy=0;};
//! destructor
    ~counter_2D(){};
//! reset the counter 
    void reset(){count=0;valuex=0;valuey=0;valuex2=0;valuey2=0;valuexy=0;};
//! add entry
    void add(
            Double_t sx, //!< x value added
            Double_t sy  //!< y value added
            ){count++; valuex+=sx; 
                       valuey+=sy; 
                       valuex2+=sx*sx;
                       valuey2+=sy*sy;
                       valuexy+=sx*sy;
             };
//! get the number of entries
    Int_t getN(){return count;};
//! get the sum of x values
    Double_t get_x(){return valuex;};
//! get the sum of x values
    Double_t get_y(){return valuey;};
//! get the sum of x2
    Double_t get_x2(){return valuex2;};
//! get the sum of y2
    Double_t get_y2(){return valuey2;};
//! get the sum of xy
    Double_t get_xy(){return valuexy;};
//! get the mean x value
    Double_t mean_x(){return valuex/count;};
//! get the mean y value
    Double_t mean_y(){return valuey/count;};
//! get the variance of x
    Double_t var_x(){return valuex2/(count-1)-valuex/count*valuex/(count-1);};
//! get the variance of y
    Double_t var_y(){return valuey2/(count-1)-valuey/count*valuey/(count-1);};
//! get the covariance
    Double_t cov(){return valuexy/(count-1)-valuex/count*valuey/(count-1);};
//! get the correlation
    Double_t corr(){return cov()/sqrt(var_x()*var_y());};

};

#endif

