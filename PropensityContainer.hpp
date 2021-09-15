/*
  Data container for internal, external and total reaction propensities.
  Total is always the sum of int and ext propensities.
  
  Int: Propensity that only depends on internal state.
  Ext: Propensity that depends on both internal and external states.
       If external state changes, all propensities change.
*/

#ifndef _PROPENSITYCONTAINER_H_
#define _PROPENSITYCONTAINER_H_

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>

#include "main.hpp"

class PropensityContainer
{
private:
    
public: //Constructors
    //Create new propensity container with length plen.
    PropensityContainer(int plen)
    {
      len = plen;
      qint = 0.;
      qA = 0;
      qB = 0;
      qKidA = 0;
      qtot = 0.;
      
      qint_list = std::vector<double>(plen, 0);
      qA_list = std::vector<double>(plen, 0);
      qB_list = std::vector<double>(plen, 0);
      qtot_list = std::vector<double>(plen, 0);
      qKidA_list = std::vector<double>(plen, 0);
    }

public: // Public function definitions.   

    int choose_index(double rnd);
   
    double sum_qint();
    double sum_qA();
    double sum_qB();
    double sum_qKidA();
    double sum_qtot();

    
    /* Setters */
    void update_qAll_el(double, double, double, double, int);
    void set_qAll_el(double, double, double, double, int);
      
    void update_qint_el(double, int);
    void update_qA_el(double, int);
    void update_qB_el(double, int);
    void update_qKidA_el(double, int);

    void set_qint_el(double, int);
    void set_qA_el(double, int);
    void set_qB_el(double, int);
    void set_qKidA_el(double, int);
    void set_qA_all(Hexamer *hexamers);
    void set_qB_all(Hexamer *hexamers);
    void set_qKidA_all(Hexamer *hexamers);
    
          
    /* Getters */
    double get_qint()
    {
      return this->qint;
    }

    double get_qA()
    {
      return this->qA;
    }
    
    double get_qtot()
    {
      return this->qtot;
    }
    
    double get_qint_el(int index)
    {
      return this->qint_list[index];
    }

    double get_qA_el(int index)
    {
      return this->qA_list[index];
    }
    
    double get_qtot_el(int index)
    {
      return this->qtot_list[index];
    }
    
    int get_len()
    {
      return this->len;
    }

private: // Private member variables.
    std::vector<double> qint_list;
    std::vector<double> qA_list;
    std::vector<double> qB_list;
    std::vector<double> qKidA_list;
    std::vector<double> qtot_list;
    
    double qint;
    double qA;
    double qB;
    double qKidA;
    double qtot;
    
    int len;
};

#endif
