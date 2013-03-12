#ifndef __COORDINATE_H__
#define __COORDINATE_H__

class Coordinate
{  
public:
  Coordinate() :row(0), col(0) {};
	
	Coordinate(double _row, double _col) :row(_row), col(_col) {};
	
	inline Coordinate operator+ (const Coordinate param) const {
  		Coordinate temp;
  		temp.row = row + param.row;
  		temp.col = col + param.col;
		return (temp);
	};
  	inline Coordinate operator+ (double param) const {
		Coordinate temp;
		temp.row = row + param;
		temp.col = col + param;
		return (temp);
	};

	inline Coordinate operator- () const {
		Coordinate temp;
		temp.row = - row;
		temp.col = - col;
		return (temp);
	};
  	inline Coordinate operator- (const Coordinate param) const {
  		Coordinate temp;
  		temp.row = row - param.row;
  		temp.col = col - param.col;
  		return (temp);
	};
	inline Coordinate operator- (double param) const {
		Coordinate temp;
		temp.row = row - param;
		temp.col = col - param;
		return (temp);
	};

	inline Coordinate operator/ (const Coordinate param) {
		Coordinate temp;
		temp.row = row / param.row;
		temp.col = col / param.col;
		return (temp);
	};
	
	inline Coordinate operator/ (double param) {
		Coordinate temp;
		temp.row = row / param;
		temp.col = col / param;
		return (temp);
	};
	
	inline Coordinate sqr (){
		Coordinate temp;
		temp.row = row * row;
		temp.col = col * col;
		return (temp);
	};
	
public:
  	double row, col;
};

#endif
