#ifndef __SMESSAGE_H__
#define __SMESSAGE_H__

#include <DTwoDimArray.h>
#include <Coordinate.h>
#include <string.h>
#include <utility>

// A very simple image class.
//
class SMessage
{
public:
  SMessage() {
		_rows = 0;
		_cols = 0;
		_childs = 0;
	}
	SMessage(int __rows, int __cols, int __childs)
	{ 
		_rows = __rows;
		_cols = __cols;
		_childs = __childs;
		
		vals = _DTwoDimArray<double>(__rows, __cols);
		coords = vector<_DTwoDimArray<Coordinate> >(__childs);
		for (int i = 0; i < __childs; i++) {
			coords[i] = _DTwoDimArray<Coordinate>(__rows, __cols);
		}
    }
    int rows() {
    	return _rows;
    }
    int cols() {
    	return _cols;
    }
public:
	int _rows, _cols, _childs;
	_DTwoDimArray<double> vals;
	vector<_DTwoDimArray<Coordinate> > coords;
};

#endif
