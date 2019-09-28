/*
 * CherenkovParameters.h
 *
 *  Created on: Feb 5, 2019
 *      Author: stenderm
 */

#ifndef SRC_CherenkovParameters_H_
#define SRC_CherenkovParameters_H_

#include <string>
#include <iostream>

#include "Tool.h"
#include "BeamStatus.h"
#include "TriggerClass.h"
#include "Detector.h"
#include "Geometry.h"
#include "Hit.h"
#include "LAPPDHit.h"

class CherenkovParameters {
public:
	CherenkovParameters(double tankRadius, double tankHeight, Position tankCenter,double refractiveIndex);
	void setParameter(int PdgCode, double particleStartEnergy, Position startVertex, Direction startDirection, Position stopVertex);
	bool createsARing();
	void numberOfPossibleHits(Position position_pmt[], double charge[], int npmt, std::map<unsigned long, std::vector<MCLAPPDHit>>* lappdhits);
	inline int numberOfPmtHit(){return _pmt_hits;}
	inline int numberOfLappdHit(){return _lappd_hits;}
	inline double returnCherenkovRadius(){return _cherenkov_radius;}
	inline Position returnStartVertex(){return _start_vertex;}
	inline double returnCherenkovAngle(){return _cherenkov_angle;}
	inline Position returnVectorAsPositionToTank(){
		Position temp(_vector_to_tank.X(), _vector_to_tank.Y(), _vector_to_tank.Z());
		return temp;}
	Position returnRingCenter();
	double distanceBetweenPositionAndRingCenter(Position pos);
	Position returnsVertexWrtTankCenter(Position vertex);
	virtual ~CherenkovParameters();

private:
	double returnMassInGeV();
	double particleVelocity();
	bool radiatesCherenkov();
	double cherenkovAngle();
	Direction vectorToTank(Position vertex, Direction direction);
	double lengthOfDirectionVector(Direction vectorToTank);
	double lengthOfPositionVector(Position vectorToTank);
	bool isWithinTank(Position vertex);
	double cherenkovRadius(Direction vector);
	double _tank_radius;
	double _tank_height;
	Position _tank_center;
	Position _start_vertex;
	Direction _start_direction;
	Position _end_vertex;
	double _refractive_index;
	double _cherenkov_angle;
	int _pdg_code;
	double _start_energy;
	Direction _vector_to_tank;
	int _pmt_hits;
	int _lappd_hits;
	double _cherenkov_radius;

};

#endif /* SRC_CherenkovParameters_H_ */
