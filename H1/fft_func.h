/*
fft_func.h

Header file for fft_func.c
 
Created by Martin Gren on 2014-10-22.
*/

#ifndef _fft_func_h
#define _fft_func_h

extern void powerspectrum(double *, double *, int);

extern void fft_freq(double *, double, int);

extern void powerspectrum_shift(double *, int);

extern void fft_freq_shift(double *, double, int);

#endif
