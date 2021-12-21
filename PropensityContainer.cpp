#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>

#include "PropensityContainer.hpp"
#include "Hexamer.hpp"


void PropensityContainer::update_qAll_el(double qint_el, double qA_el, double qB_el, double qKidA_el, int index)
{
  qint -= qint_list[index];
  qint += qint_el;
  
  qA -= qA_list[index];
  qA += qA_el;

  qB -= qB_list[index];
  qB += qB_el;

  qKidA -= qKidA_list[index];
  qKidA += qKidA_el;

  qtot = qint + qA + qB + qKidA;
 
  qint_list[index] = qint_el;
  qA_list[index] = qA_el;
  qB_list[index] = qB_el;
  qKidA_list[index] = qKidA_el;
  qtot_list[index] = qint_list[index] + qA_list[index] + qB_list[index] + qKidA_list[index];

  if(qint <= 0) {
  
    sum_qint();
    sum_qtot();
  }
  if(qA <= 0) {

    sum_qA();
    sum_qtot();
  }
  if(qB <= 0) {

    sum_qB();
    sum_qtot();
  }
  if(qKidA < 0) {

    sum_qKidA();
    sum_qtot();
  }
}


void PropensityContainer::set_qAll_el(double qint_el, double qA_el, double qB_el, double qKidA_el, int index)
{ 
  qint_list[index] = qint_el;
  sum_qint();
  
  qA_list[index] = qA_el;
  sum_qA();
  
  qB_list[index] = qB_el;
  sum_qB();

  qKidA_list[index] = qKidA_el;
  sum_qKidA();

  qtot_list[index] = qint_list[index] + qA_list[index] + qB_list[index] + qKidA_list[index];
  sum_qtot();
}


void PropensityContainer::update_qint_el(double qint_el, int index)
{
  qint -= qint_list[index];
  qint += qint_el;
  qtot = qint + qA + qB + qKidA;
    
  qint_list[index] = qint_el;
  qtot_list[index] = qint_list[index] + qA_list[index] + qB_list[index] + qKidA_list[index];
  
  if(qint < 0) {

    sum_qint();
    sum_qtot();
  }
}


void PropensityContainer::set_qint_el(double qint_el, int index)
{      
  qint_list[index] = qint_el;
  sum_qint();
 
  qtot_list[index] = qint_list[index] + qA_list[index] + qB_list[index];
  sum_qtot();
}
  
    
void PropensityContainer::update_qA_el(double qA_el, int index)
{     
  qA -= qA_list[index];
  qA += qA_el;
  qtot = qint + qA + qB + qKidA;
    
  qA_list[index] = qA_el;
  qtot_list[index] = qint_list[index] + qA_list[index] + qB_list[index] + qKidA_list[index];
}


void PropensityContainer::set_qA_el(double qA_el, int index)
{     
  qA_list[index] = qA_el;
  sum_qA();
 
  qtot_list[index] = qint_list[index] + qA_list[index] + qB_list[index] + qKidA_list[index];
  sum_qtot();
}

void PropensityContainer::update_qB_el(double qB_el, int index)
{     
  qB -= qB_list[index];
  qB += qB_el;
  qtot = qint + qA + qB + qKidA;
    
  qB_list[index] = qB_el;
  qtot_list[index] = qint_list[index] + qA_list[index] + qB_list[index] + qKidA_list[index];
}

void PropensityContainer::set_qB_el(double qB_el, int index)
{     
  qB_list[index] = qB_el;
  sum_qB();
 
  qtot_list[index] = qint_list[index] + qA_list[index] + qB_list[index] + qKidA_list[index];
  sum_qtot();     
}

void PropensityContainer::update_qKidA_el(double qKidA_el, int index)
{
	qKidA -= qKidA_list[index];
	qKidA += qKidA_el;
	qtot = qint + qA + qB + qKidA;
	
	qKidA_list[index] = qKidA_el;
	qtot_list[index] = qint_list[index] + qA_list[index] + qB_list[index] + qKidA_list[index];
}

void PropensityContainer::set_qKidA_el(double qKidA_el, int index)
{
	qKidA_list[index] = qKidA_el;
	sum_qKidA();

	qtot_list[index] = qint_list[index] + qA_list[index] + qB_list[index] + qKidA_list[index];
	sum_qtot();     
}

void PropensityContainer::set_qA_all(Hexamer *hexamers)
{     
  for(int i(0); i<len; i++)
  {
    qA_list[i] = hexamers[i].update_KaiA_prop();
    qtot_list[i] = qint_list[i] + qA_list[i] + qB_list[i] + qKidA_list[i];
  }
      
  sum_qA();
  sum_qtot();       
}

void PropensityContainer::set_qB_all(Hexamer *hexamers)
{     
  for(int i(0); i<len; i++)
  {
    qB_list[i] = hexamers[i].update_KaiB_prop();
    qtot_list[i] = qint_list[i] + qA_list[i] + qB_list[i] + qKidA_list[i];
  }
      
  sum_qB();
  sum_qtot();       
}

void PropensityContainer::set_qKidA_all(Hexamer *hexamers)
{
	for(int i = 0; i < len; i += 1) {

		qKidA_list[i] = hexamers[i].update_KidA_prop();
		qtot_list[i] = qint_list[i] + qA_list[i] + qB_list[i] + qKidA_list[i];
	}

	sum_qKidA();
	sum_qtot();
}

double PropensityContainer::sum_qint()
{
  double sum(0.0);
    
  for(int i(0); i<len; i++)
  {
    sum += qint_list[i];
  }
      
  return this->qint = sum;
}
    
    
double PropensityContainer::sum_qA()
{
  double sum(0.0);
    
  for(int i(0); i<len; i++)
  {
    sum += qA_list[i];
  }
      
  return this->qA = sum;
}


double PropensityContainer::sum_qB()
{
  double sum(0.0);
    
  for(int i(0); i<len; i++)
  {
    sum += qB_list[i];
  }
      
  return this->qB = sum;
}


double PropensityContainer::sum_qtot()
{
  double sum(0.0);
    
  for(int i(0); i<len; i++)
  {
    sum += qtot_list[i];
  }
      
  return this->qtot = sum;
}

double PropensityContainer::sum_qKidA()
{

	double sum(0);

	for(int i = 0; i < len; i += 1) {

		sum += qKidA_list[i];
	}
	
	this->qKidA = sum;
	return sum;
}

int PropensityContainer::choose_index(double rnd)
{
  double qacc(qtot_list[0]);
  int i(0);
  
  rnd*=SCALE;  
  rnd*=qtot;
  while( rnd > qacc && i < len - 1)
  {
    i += 1;
    qacc += qtot_list[i];
  }
  
  return i;
}
