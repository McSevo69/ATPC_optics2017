#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "Maxfiles.h"
#include <MaxSLiCInterface.h>

typedef float dataType;

struct point {
	float x;
	float y;
	float z;
	float x4;
	float x5;
	float x6;
	float x7;
	float x8;
	float cd; //core distance
	float rd; //reachability distance
	int processed;
};

int checkResults(size_t *order1, size_t *order2) {
	for (size_t i=0; i<Vectors_maxN; i++) if (order1[i] != order2[i]) return 1;

	return 0;
}

void printInData(struct point *points) {
	printf("\n");
	for (size_t i = 0; i < Vectors_maxN; i++) {		
		printf("%f %f %f %f %f %f %f %f", points[i].x, points[i].y, points[i].z,
			points[i].x4, points[i].x5, points[i].x6, points[i].x7, points[i].x8);              
		printf("\n");
	}
}

//adapted from http://c-faq.com/lib/gaussian.html
float gaussrand() {
	static float V1, V2, S;
	static int phase = 0;
	float X;

	if(phase == 0) {
		do {
			float U1 = (float)rand() / RAND_MAX;
			float U2 = (float)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}


void createRandomData(struct point *dataset, int size) {
                
        float meanX, meanY, meanZ, meanX4, meanX5, meanX6, meanX7, meanX8; 
        float variance;      
                
        for (int i = 0; i < size; i++) {
		meanX = 30 + (120 - 30) * (float)rand() / RAND_MAX;
		meanY = 30 + (120 - 30) * (float)rand() / RAND_MAX;
		meanZ = 30 + (120 - 30) * (float)rand() / RAND_MAX;
		meanX4 = 30 + (120 - 30) * (float)rand() / RAND_MAX;
		meanX5 = 30 + (120 - 30) * (float)rand() / RAND_MAX;
		meanX6 = 30 + (120 - 30) * (float)rand() / RAND_MAX;
		meanX7 = 30 + (120 - 30) * (float)rand() / RAND_MAX;
		meanX8 = 30 + (120 - 30) * (float)rand() / RAND_MAX;
		variance = 6 + (10 - 6) * (float)rand() / RAND_MAX;

		dataset[i].x = gaussrand() * variance + meanX;
		dataset[i].y = gaussrand() * variance + meanY;
		dataset[i].z = gaussrand() * variance + meanZ;
		dataset[i].x4 = gaussrand() * variance + meanX4;
		dataset[i].x5 = gaussrand() * variance + meanX5;
		dataset[i].x6 = gaussrand() * variance + meanX6;
		dataset[i].x7 = gaussrand() * variance + meanX7;
		dataset[i].x8 = gaussrand() * variance + meanX8;
		dataset[i].rd = INFINITY;
		dataset[i].cd = INFINITY;
		dataset[i].processed = 0;
        }
}

struct point* loadDataFromFile(char *filename) {

	char *line = NULL;
	size_t n = 0, i = 0;

	struct point *inPointsBuf = malloc(Vectors_maxN * sizeof *inPointsBuf);	
	
	FILE *dataFile = fopen(filename, "r");

	while(getline(&line, &n, dataFile) != -1) {
		int items = sscanf(line, "%f %f %f %f %f %f %f %f", &inPointsBuf[i].x, &inPointsBuf[i].y, &inPointsBuf[i].z,
		&inPointsBuf[i].x4, &inPointsBuf[i].x5, &inPointsBuf[i].x6, &inPointsBuf[i].x7, &inPointsBuf[i].x8);
		inPointsBuf[i].rd = INFINITY;
		inPointsBuf[i].cd = INFINITY;
		inPointsBuf[i].processed = 0;

		if (items != 8) {
			printf("Invalid input file\n");
			exit(1);
		}
		i++;
	}

	fclose(dataFile);

	return inPointsBuf;
}

__inline__ float getEuclidianDistance(struct point point1, struct point point2) {
	return sqrt(pow(point1.x-point2.x, 2) + pow(point1.y-point2.y, 2) + pow(point1.z-point2.z, 2)
		+ pow(point1.x4-point2.x4, 2) + pow(point1.x5-point2.x5, 2) + pow(point1.x6-point2.x6, 2)
		+ pow(point1.x7-point2.x7, 2) + pow(point1.x8-point2.x8, 2));
}

__inline__ float getEuclidianDistanceLUT(size_t p1, size_t p2, dataType *LUT) {
	return (p1 > p2) ? LUT[p2*Vectors_maxN+p1] : LUT[p1*Vectors_maxN+p2];
}


void sortRDAscending(struct point *list, int length) {

	for (int i = 0; i < length; i++) {
		for (int j = 0; j < length; j++) {
			if (list[j].rd > list[i].rd) {
				struct point temp = list[i];
				list[i] = list[j];
				list[j] = temp;
			}
		}
	}
}

float getCoreDistance(struct point *dataSet, struct point p, int minPts, float eps) {

	struct point *distList = malloc(Vectors_maxN*sizeof *distList);

	for (size_t i=0; i<Vectors_maxN; i++) {	
		distList[i].rd = getEuclidianDistance(p, dataSet[i]);	
	}

	sortRDAscending(distList, Vectors_maxN);

	float cd = (distList[minPts].rd <= eps) ? distList[minPts].rd : INFINITY; 

	free(distList);
	return cd;
}


float getCoreDistanceLUT(size_t pointIndex, dataType *euclidianLUT, int minPts, float eps) {

	struct point *distList = malloc(Vectors_maxN*sizeof *distList);

	for (size_t i=0; i<Vectors_maxN; i++) {	
		distList[i].rd = getEuclidianDistanceLUT(pointIndex, i, euclidianLUT);
	}

	sortRDAscending(distList, Vectors_maxN);

	float cd = (distList[minPts].rd <= eps) ? distList[minPts].rd : INFINITY; 

	free(distList);
	return cd;
}


int getNeighbors(struct point *dataSet, size_t *neighborIndices, size_t p, float eps) {

	int cnt = 0;

	for (size_t i=0; i<Vectors_maxN; i++) {
		if (i != p && getEuclidianDistance(dataSet[p], dataSet[i]) <= eps) {
			neighborIndices[cnt++] = i;
		}
	}
	
	return cnt;
}

int getNeighborsLUT(size_t *neighborIndices, dataType *euclidianLUT, size_t p, float eps) {

	int cnt = 0;

	for (size_t i=0; i<Vectors_maxN; i++) {
		if (i != p && getEuclidianDistanceLUT(p, i, euclidianLUT) <= eps) {
			neighborIndices[cnt++] = i;
		}
	}
	
	return cnt;
}


//ordered by rd
void insertIntoList(struct point *dataSet, size_t *list, size_t p, int n) {
	
	if (n == 0) list[0] = p; //first element
	else if (dataSet[p].rd > dataSet[list[n-1]].rd) list[n] = p; //worst rd
	else {
		
		int i=-1;
		while (dataSet[p].rd > dataSet[list[++i]].rd); //new place found		

		for (int j=n; j >= i; j--) { //moving elements one place backwards
			list[j+1] = list[j];	
		}

		list[i] = p;
	}
}

void removeFromList(size_t *list, size_t p, int n) {
	if (list[n-1] != p) {
		int i=0;
		while (list[i++] != p); //

		for (i; i < n; i++) { //moving elements one place forward
			list[i-1] = list[i];	
		}
	}
}

int moveUp(struct point *dataSet, size_t *list, size_t p, int n) {
	int isMember = 0;

	for (int i=0; i<n; i++) {
		if (list[i] == p) {
			isMember = 1;
			break;
		} 		
	}

	if (isMember) {
		removeFromList(list, p, n--);
		insertIntoList(dataSet, list, p, n++);
		return n;
	} else {
		insertIntoList(dataSet, list, p, n++);
		return n;
	}
}

void nextInQueue(size_t *queue, int n) {	
	for (int i=0; i<n-1; i++) queue[i] = queue[i+1];
	//queue[n-1] = (int*) NULL;
}


int update(struct point* dataSet, size_t *seedList, size_t *neighbors, size_t p, int neighborCnt, int cnt) {
	float coredist = dataSet[p].cd;
	
	for (int i=0; i<neighborCnt; i++) {
		if (!dataSet[neighbors[i]].processed) {
			float newRD = fmax(coredist, getEuclidianDistance(dataSet[p], dataSet[neighbors[i]]));
			if (isinf(dataSet[neighbors[i]].rd)) { //-> no member of seedlist
				dataSet[neighbors[i]].rd = newRD;
				insertIntoList(dataSet, seedList, neighbors[i], cnt++);
			} else {
				if (newRD < dataSet[neighbors[i]].rd) {
					dataSet[neighbors[i]].rd = newRD;
					cnt = moveUp(dataSet, seedList, neighbors[i], cnt);
				}
			}
		}
	}
	return cnt;	
}

int updateLUT(struct point* dataSet, size_t *seedList, size_t *neighbors, dataType *euclidianLUT, size_t p, int neighborCnt, int cnt) {
	float coredist = dataSet[p].cd;
	
	for (int i=0; i<neighborCnt; i++) {
		if (!dataSet[neighbors[i]].processed) {
			float newRD = fmax(coredist, getEuclidianDistanceLUT(p, neighbors[i], euclidianLUT));
			if (isinf(dataSet[neighbors[i]].rd)) { //-> no member of seedlist
				dataSet[neighbors[i]].rd = newRD;
				insertIntoList(dataSet, seedList, neighbors[i], cnt++);
			} else {
				if (newRD < dataSet[neighbors[i]].rd) {
					dataSet[neighbors[i]].rd = newRD;
					cnt = moveUp(dataSet, seedList, neighbors[i], cnt);
				}
			}
		}
	}
	return cnt;	
}

//main-algorithm
void doOpticsCPU(struct point *dataSet, size_t *output, float eps, int minPts) {

	int cnt = 0, queueCnt = 0, neighborCnt = 0;
	//both for storing indices
	size_t *queue = malloc(Vectors_maxN * sizeof(size_t));
	size_t *neighbors = malloc(Vectors_maxN * sizeof(size_t)); 

	for (size_t i=0; i<Vectors_maxN; i++) {
		if (!dataSet[i].processed) {			
			neighborCnt = getNeighbors(dataSet, neighbors, i, eps);
			dataSet[i].cd = getCoreDistance(dataSet, dataSet[i], minPts, eps);
			dataSet[i].processed = 1;
			output[cnt++] = i;

			if (isfinite(dataSet[i].cd)) {
				queueCnt = update(dataSet, queue, neighbors, i, neighborCnt, queueCnt);

				while (queueCnt > 0) {
					size_t next = queue[0];
					nextInQueue(queue, queueCnt--); //removing first element

					neighborCnt = getNeighbors(dataSet, neighbors, next, eps);
					dataSet[next].processed = 1;
					output[cnt++] = next;
					
					dataSet[next].cd = getCoreDistance(dataSet, dataSet[next], minPts, eps);
					if (isfinite(dataSet[next].cd)) {
						update(dataSet, queue, neighbors, next, neighborCnt, queueCnt);
					}					
				}
			}			
		}
	}

	free(queue);
	free(neighbors);
}


//main-algorithm
void doOpticsLUT(struct point *dataSet, size_t *output, dataType *euclidianLUT, float eps, int minPts) {

	int cnt = 0, queueCnt = 0, neighborCnt = 0;
	//both for storing indices
	size_t *queue = malloc(Vectors_maxN * sizeof(size_t));
	size_t *neighbors = malloc(Vectors_maxN * sizeof(size_t)); 

	for (size_t i=0; i<Vectors_maxN; i++) {
		if (!dataSet[i].processed) {			
			neighborCnt = getNeighborsLUT(neighbors, euclidianLUT, i, eps);
			dataSet[i].cd = getCoreDistanceLUT(i, euclidianLUT, minPts, eps);
			dataSet[i].processed = 1;
			output[cnt++] = i;

			if (isfinite(dataSet[i].cd)) {
				queueCnt = updateLUT(dataSet, queue, neighbors, euclidianLUT, i, neighborCnt, queueCnt);

				while (queueCnt > 0) {
					size_t next = queue[0];
					nextInQueue(queue, queueCnt--); //removing first element

					neighborCnt = getNeighborsLUT(neighbors, euclidianLUT, next, eps);
					dataSet[next].processed = 1;
					output[cnt++] = next;
					
					dataSet[next].cd = getCoreDistanceLUT(next, euclidianLUT, minPts, eps);
					if (isfinite(dataSet[next].cd)) {
						updateLUT(dataSet, queue, neighbors, euclidianLUT, next, neighborCnt, queueCnt);
					}					
				}
			}			
		}
	}

	free(queue);
	free(neighbors);
}

void exportReachability(struct point *dataSet, size_t *outputOrder, char *filename) {

	FILE* results;
	
	results = fopen(filename, "w");

	for (size_t i=0; i<Vectors_maxN; i++) {
		fprintf(results, "%f %f %f %f %f %f %f %f %f\n", dataSet[outputOrder[i]].x, dataSet[outputOrder[i]].y,
			dataSet[outputOrder[i]].z, dataSet[outputOrder[i]].x4, dataSet[outputOrder[i]].x5,
			dataSet[outputOrder[i]].x6, dataSet[outputOrder[i]].x7, dataSet[outputOrder[i]].x8,
			dataSet[outputOrder[i]].rd);
	}
	
	fclose(results);
}


int convertArgToInt(char * str) {
	if (strcmp ("-in", str) == 0) return 1;
	else if (strcmp ("-out", str) == 0) return 2;
	else if (strcmp ("-print", str) == 0) return 3;
	else if (strcmp ("-benchmark", str) == 0) return 4;
	else if (strcmp ("-minPts", str) == 0) return 5;
	else if (strcmp ("-eps", str) == 0) return 6;
	else return 0;
}

int main(int argc, char *argv[]) {
	//time measurement
	struct timeval begin, end;
	double timeSpentCPU = 0, timeSpent = 0;

	struct point *inPoints, *inPointsDFE, *outPoints;
	size_t *outputOrder, *outputOrderDFE;
	
	printf("=================\n");
	printf("OPTICS stream\n");
	printf("Precision: %ld\n",sizeof(dataType)*8);
	printf("=================\n");
	
	size_t benchmark = 0, print = 0, fileGiven = 0;
	int minPts = 3; 
	float eps = 20.0;
	char *inPath, *outPath = "rd.csv";

	//commandline parameter parsing
	for (int i=1; i<argc; i++) {
		switch(convertArgToInt(argv[i])) {
			case 1: inPath = argv[++i]; fileGiven = 1; break;
			case 2: outPath = argv[++i]; break;
			case 3: print = 1; break;
			case 4: benchmark = 1; break;
			case 5: minPts = atoi(argv[++i]); break;
			case 6: eps = atof(argv[++i]); break;
			default: break;
		}
	}

	if (fileGiven) {
		printf("Loading data from file %s ... \n", inPath);
		inPointsDFE = loadDataFromFile(inPath);
		inPoints = malloc(Vectors_maxN * sizeof *inPoints);
		for (size_t i=0; i<Vectors_maxN; i++) inPoints[i] = inPointsDFE[i];
	} else {
		printf("Generating input data.\n");
		inPointsDFE = malloc(Vectors_maxN * sizeof *inPoints);
		createRandomData(inPointsDFE, Vectors_maxN);
		inPoints = malloc(Vectors_maxN * sizeof *inPoints);
		for (size_t i=0; i<Vectors_maxN; i++) inPoints[i] = inPointsDFE[i];
	}

	if (print) {
		printf("Printing input data.\n");
		printInData(inPoints);
	}

	outPoints = malloc(Vectors_maxN * sizeof *outPoints);
	dataType *euclidianLUT = malloc(Vectors_maxN * Vectors_maxN * sizeof(dataType));

	outputOrder = malloc(Vectors_maxN * sizeof(size_t));
	outputOrderDFE = malloc(Vectors_maxN * sizeof(size_t));	

	if (benchmark) {
		outputOrder = malloc(Vectors_maxN * sizeof(size_t));
		gettimeofday(&begin, NULL);
		doOpticsCPU(inPoints, outputOrder, eps, minPts);
		gettimeofday(&end, NULL);
		timeSpentCPU += (end.tv_sec - begin.tv_sec) +
		    ((end.tv_usec - begin.tv_usec)/1000000.0);
		printf("Time CPU: %f\n", timeSpentCPU);
	}	

	dataType *inValuesDFE = calloc(8*Vectors_maxN, sizeof(dataType));
	
	for (size_t i=0; i<Vectors_maxN; i++) {
		inValuesDFE[i] = inPointsDFE[i].x;
		inValuesDFE[Vectors_maxN+i] = inPointsDFE[i].y;
		inValuesDFE[2*Vectors_maxN+i] = inPointsDFE[i].z;
		inValuesDFE[3*Vectors_maxN+i] = inPointsDFE[i].x4;
		inValuesDFE[4*Vectors_maxN+i] = inPointsDFE[i].x5;
		inValuesDFE[5*Vectors_maxN+i] = inPointsDFE[i].x6;
		inValuesDFE[6*Vectors_maxN+i] = inPointsDFE[i].x7;
		inValuesDFE[7*Vectors_maxN+i] = inPointsDFE[i].x8;
	}

	gettimeofday(&begin, NULL);
	Vectors(Vectors_maxN, inValuesDFE, &inValuesDFE[Vectors_maxN], &inValuesDFE[2*Vectors_maxN], 
		&inValuesDFE[3*Vectors_maxN], &inValuesDFE[4*Vectors_maxN], &inValuesDFE[5*Vectors_maxN],
		&inValuesDFE[6*Vectors_maxN], &inValuesDFE[7*Vectors_maxN], euclidianLUT);
	doOpticsLUT(inPointsDFE, outputOrderDFE, euclidianLUT, eps, minPts);
	gettimeofday(&end, NULL);
	timeSpent += (end.tv_sec - begin.tv_sec) +
            ((end.tv_usec - begin.tv_usec)/1000000.0);

	printf("Time DFE: %f\n", timeSpent);

	//Exporting benchmark results
	if (benchmark) {
		int result = checkResults(outputOrderDFE, outputOrder);
		if (result) printf("Test failed.\n");
		else printf("Test succeeded.\n");

		FILE* time_res;
		char filename[32];
		printf("Exporting to file...\n");
		snprintf(filename, sizeof(filename), "benchmark.csv");
		time_res = fopen(filename,"a");
		fprintf(time_res, "%d, %f, %f\n", Vectors_maxN, timeSpentCPU, timeSpent);
		fclose(time_res);
	}
	
	printf("Exporting results...\n");
	exportReachability(inPoints, outputOrderDFE, outPath);

	free(inPoints);
	free(inPointsDFE);
	free(inValuesDFE);
	free(outPoints);
	free(outputOrder);
	free(outputOrderDFE);
	
	return 0;
}
