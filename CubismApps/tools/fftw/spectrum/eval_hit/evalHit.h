/* File:   evalHit.h */
/* Date:   July 2015*/
/* Author: Ursula Rasthofer */
/* Tag:    evaluations for hit */
/* Copyright 2015 ETH Zurich. All Rights Reserved. */
#pragma once

#include "../evalBase.h"
#include "../../FFTWVectorMPI.h"
#include <Cubism/ArgumentParser.h>
using namespace cubism;

class evalHit: public evalBase, public FFTWVectorMPI
{

public:

  // constructor
  evalHit(ArgumentParser& myparser);
  // destructor
  ~evalHit() {}

  // compute initial field
  void compute();

  // clean up and free memory
  void dispose();

private:

  // helpers
  ArgumentParser& myparser_;

  // input structure
  Real* data_2;

};
