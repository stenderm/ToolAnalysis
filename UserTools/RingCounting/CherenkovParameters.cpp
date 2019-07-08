/*
 * CherenkovParameters.cpp
 *
 *  Created on: Feb 5, 2019
 *      Author: stenderm
 */

#include "CherenkovParameters.h"

CherenkovParameters::CherenkovParameters(double tankRadius, double tankHeight,
		Position tankCenter, double refractiveIndex) {
	_tank_center = tankCenter;
	_tank_radius = tankRadius;
	_tank_height = tankHeight;
	_refractive_index = refractiveIndex;
}

CherenkovParameters::~CherenkovParameters() {
	// TODO Auto-generated destructor stub
}

void CherenkovParameters::setParameter(int PdgCode, double particleStartEnergy,
		Position startVertex, Direction startDirection, Position stopVertex) {
	_pdg_code = PdgCode;
	_start_energy = particleStartEnergy;
	_start_vertex = returnsVertexWrtTankCenter(startVertex);
	_start_direction = startDirection;
	_end_vertex = returnsVertexWrtTankCenter(stopVertex);
	if (createsARing()) {
		_cherenkov_angle = cherenkovAngle();
		_vector_to_tank = vectorToTank(_start_vertex, _start_direction);
		_cherenkov_radius = cherenkovRadius(_vector_to_tank);
	}
}

bool CherenkovParameters::createsARing() {
	if (radiatesCherenkov()) {
		if (isWithinTank(_start_vertex) || isWithinTank(_end_vertex)) {
			return true;
		} else {
			return false;
		}
	} else {
		return false;
	}
}

//Returns mass of particles. Returns zero for an invalid Code and for Photons, Gammas, Neutrinos and for uncharged particles, as they do not do Cherenkov radiation
double CherenkovParameters::returnMassInGeV() {
	switch (abs(_pdg_code)) {
	case 11:
		return 0.000511;
	case 13:
		return 0.105658;
	case 15:
		return 1.77686;
	case 211:
		return 0.139570;
	case 321:
		return 0.493677;
	case 2212:
		return 0.938272;
	case 3122:
		return 1.115638;
	case 3112:
		return 1.197449;
	case 3222:
		return 1.18937;
	case 3312:
		return 1.32171;
	case 3334:
		return 1.67245;
	case 3328:
		return 0.372737;
	default:
		return 0;
	}
}

double CherenkovParameters::particleVelocity() {
	double velocity = TMath::Sqrt(
			1 - pow(returnMassInGeV() / _start_energy, 2));
	return velocity;
}

bool CherenkovParameters::radiatesCherenkov() {
	double velocity = particleVelocity();
	if (velocity > (1 / _refractive_index) && returnMassInGeV() != 0) {
		return true;
	} else {
		return false;
	}
}

Position CherenkovParameters::returnsVertexWrtTankCenter(Position vertex) {
	Position vertexWrtToTankCenter;
	vertexWrtToTankCenter.SetX((vertex.X() - _tank_center.X()));
	vertexWrtToTankCenter.SetY((vertex.Y() - _tank_center.Y()));
	vertexWrtToTankCenter.SetZ((vertex.Z() - _tank_center.Z()));
	return vertexWrtToTankCenter;
}

double CherenkovParameters::cherenkovAngle() {
	double angle = acos(1 / (_refractive_index * particleVelocity()));
	return angle;
}

Direction CherenkovParameters::vectorToTank(Position vertex,
		Direction directionVector) {
	double directionX = directionVector.X();
	double directionY = directionVector.Y();
	double directionZ = directionVector.Z();
	double realVertexX = vertex.X();
	double realVertexY = vertex.Y();
	double realVertexZ = vertex.Z();
	double denominator = ((directionX * directionX) + (directionZ * directionZ));
	double p = 2
			* (((realVertexX * directionX) + (realVertexZ * directionZ))
					/ denominator);
	double q = (((realVertexZ * realVertexZ) + (realVertexX * realVertexX)
			- (_tank_radius * _tank_radius)) / denominator);
	double lambda1 = ((-p / 2) + TMath::Sqrt((p * p / 4) - q));
	double lambda2 = ((-p / 2) - TMath::Sqrt((p * p / 4) - q));
	double lambdaRadius = 0;
	if (lambda1 > 0) {
		lambdaRadius = lambda1;
	} else {
		lambdaRadius = lambda2;
	}
	double lambdaHeight = 0.0;
	lambdaHeight = (_tank_height - realVertexY) / directionY;
	// p = 2*realVertexY/directionY;
	// q = ((realVertexY*realVertexY)-(_tank_height*_tank_height))/(directionY*directionY);
	// lambda1 = (-(1/2)*p)+sqrt((p*p/4)-q);
	// lambda2 = (-(1/2)*p)-sqrt((p*p/4)-q);
	// double lambdaHeight = 0.0;
	// if(lambda1 > 0){
	//   lambdaHeight = lambda1;
	// }
	// else{
	//   lambdaHeight = lambda2;
	// }
	double lambdaRight = 0;
	if (abs(lambdaHeight) > lambdaRadius) {
		lambdaRight = lambdaRadius;
	} else {
		lambdaRight = lambdaHeight;
	}
	Direction direction;
	direction.SetX((lambdaRight * directionX)); //vorher mit realVertexX+lambdaRight*directionX
	direction.SetY((lambdaRight * directionY));
	direction.SetZ((lambdaRight * directionZ));
	return direction;
}

double CherenkovParameters::lengthOfDirectionVector(Direction vector) {
	return TMath::Sqrt(
			pow(vector.X(), 2) + pow(vector.Y(), 2) + pow(vector.Z(), 2));
}

double CherenkovParameters::lengthOfPositionVector(Position vector) {
	return TMath::Sqrt(
			pow(vector.X(), 2) + pow(vector.Y(), 2) + pow(vector.Z(), 2));
}

bool CherenkovParameters::isWithinTank(Position vertex) {
	double distanceToCentralAxis = sqrt(
			vertex.X() * vertex.X() + vertex.Z() * vertex.Z());
	if (distanceToCentralAxis > _tank_radius
			|| abs(vertex.Y()) > _tank_height) {
		return false;
	} else {
		return true;
	}
}

double CherenkovParameters::cherenkovRadius(Direction vector) {
	double radius = sin(_cherenkov_angle) * lengthOfDirectionVector(vector)
			/ sin((M_PI / 2) - _cherenkov_angle);
	return radius;
}

void CherenkovParameters::numberOfPossibleHits(Position position_pmt[],
		double charge[], int npmt,
		std::map<unsigned long, std::vector<MCLAPPDHit>>* lappdhits) {
	int counterPMT = 0;
	int counterLAPPD = 0;
	for (int i_pmt = 0; i_pmt < 200; i_pmt++) {
		if (charge[i_pmt] > 1) {
			Direction vectorFromVertexToPmt(
					position_pmt[i_pmt].X() - _start_vertex.X(),
					position_pmt[i_pmt].Y() - _start_vertex.Y(),
					position_pmt[i_pmt].Z() - _start_vertex.Z());
			double dotProduct = vectorFromVertexToPmt.X() * _vector_to_tank.X()
					+ vectorFromVertexToPmt.Y() * _vector_to_tank.Y()
					+ vectorFromVertexToPmt.Z() * _vector_to_tank.Z();
			double angle = acos(
					dotProduct
							/ (lengthOfDirectionVector(_vector_to_tank)
									* lengthOfDirectionVector(
											vectorFromVertexToPmt)));
			if (angle <= _cherenkov_angle) {
				// std::cout << "pmt (" << xpmt[i_pmt] << "|" << ypmt[i_pmt] << "|" << zpmt[i_pmt] << ")" << std::endl;
				// std::cout << "charge " << charge[i_pmt] << std::endl;
				// std::cout << "Distance between PMT and Ring Center " << distanceBetweenPositionAndRingCenter(pmtPos) << std::endl;
				counterPMT++;
			}
		}
	}
	if (lappdhits) {
		for (auto&& achannel : *lappdhits) {
			unsigned long chankey = achannel.first;
			auto& hits = achannel.second;
			for (auto&& ahit : hits) {
				std::vector<double> temp_pos = ahit.GetPosition();
				Position lappdPos(temp_pos.at(0), temp_pos.at(1),
						temp_pos.at(2));
				Position lappdPosWrtCenter = returnsVertexWrtTankCenter(
						lappdPos);
				Direction vectorFromVertexToLappd;
				vectorFromVertexToLappd.SetX(
						lappdPosWrtCenter.X() - _start_vertex.X());
				vectorFromVertexToLappd.SetY(
						lappdPosWrtCenter.Y() - _start_vertex.Y());
				vectorFromVertexToLappd.SetZ(
						lappdPosWrtCenter.Z() - _start_vertex.Z());
				double dotProduct = vectorFromVertexToLappd.X()
						* _vector_to_tank.X()
						+ vectorFromVertexToLappd.Y() * _vector_to_tank.Y()
						+ vectorFromVertexToLappd.Z() * _vector_to_tank.Z();
				double angle = acos(
						dotProduct
								/ (lengthOfDirectionVector(_vector_to_tank)
										* lengthOfDirectionVector(
												vectorFromVertexToLappd)));
				if (angle <= _cherenkov_angle) {
					counterLAPPD++;
				}
			}
		}
	}
	_pmt_hits = counterPMT;
	// cout << "pmthits " << _pmt_hits << endl;
	_lappd_hits = counterLAPPD;
}

Position CherenkovParameters::returnRingCenter() {
	Position tempPos;
	tempPos.SetX(_start_vertex.X() + _vector_to_tank.X());
	tempPos.SetY(_start_vertex.Y() + _vector_to_tank.Y());
	tempPos.SetZ(_start_vertex.Z() + _vector_to_tank.Z());
	return tempPos;
}

double CherenkovParameters::distanceBetweenPositionAndRingCenter(Position pos) {
	Position tempPos;
	Position center = returnRingCenter();
	tempPos.SetX(pos.X() - center.X());
	tempPos.SetY(pos.Y() - center.Y());
	tempPos.SetZ(pos.Z() - center.Z());
	return lengthOfPositionVector(tempPos);
}
