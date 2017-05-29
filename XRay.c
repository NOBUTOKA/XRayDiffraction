#include <stdio.h>
#include <math.h>

#define PI 3.1415
#define LAMBDA 0.154

enum{
	NACL = 0,
	ZNS,
	DIAMOND,
	BCC,
	STRUCT_NUM	// Number of considered crystal structure.
};

double printPeakToHKL(int* peaks, double* HKL, int peaksNum, int HKLNum);

int main(void) {
	printf("Program Started.\n");
	// List of Sqrt(h^2+k^2+l^2) that can be observed peak.
	double HKL[STRUCT_NUM][36] = {
		{sqrt(3), sqrt(4), sqrt(8), sqrt(11), sqrt(12),
		 sqrt(16), sqrt(19), sqrt(20), sqrt(24), sqrt(27), sqrt(32), sqrt(35), sqrt(36)},
		{sqrt(3), sqrt(4), sqrt(8), sqrt(11), sqrt(12),
		 sqrt(16), sqrt(19), sqrt(20), sqrt(24), sqrt(27), sqrt(32), sqrt(35), sqrt(36)},
		{sqrt(3), sqrt(8), sqrt(11), sqrt(12), sqrt(16), sqrt(17), 
		 sqrt(19), sqrt(24), sqrt(27), sqrt(32), sqrt(35)},
		{sqrt(2), sqrt(4), sqrt(6), sqrt(8), sqrt(10), sqrt(12), sqrt(14), sqrt(16),
		 sqrt(18), sqrt(20), sqrt(22), sqrt(24), sqrt(26), sqrt(30), sqrt(32), sqrt(34), sqrt(36)}};
	int HKLNum[STRUCT_NUM] = {13, 13, 11, 17}; // Size of CSF[i];
	int peaksNum;	// Number of actual peaks.
	scanf("%d\n", &peaksNum);
	int peaks[peaksNum];	// Double of peaks degree(2*theta).
	for(int i = 0; i < peaksNum; i++){
		scanf("%d\n", &peaks[i]);
	}

	/* printPeakToHKL(peaks, HKL[NACL], peaksNum, HKLNum[NACL]);  */
	/* return 0; */
	
	double maxInvStddev = 0;
	int mostPossibleStruct[2] = {};
	double stddevInMostPossible = 0;
	for(int j = 0; j < 10; j++){
		for(int i = 0; i < STRUCT_NUM; i++){
			switch(i){
			case NACL:
				printf("NaCl:\n");
				break;
			case ZNS:
				printf("ZnS:\n");
				break;
			case DIAMOND:
				printf("Diamond:\n");
				break;
			case BCC:
				printf("Bcc:\n");
				break;
			}
		
			double stddev = printPeakToHKL(peaks, HKL[i] + j, peaksNum, HKLNum[i] - j);
			if(maxInvStddev < (1.0 / stddev)){
				maxInvStddev = 1.0 / stddev;
				stddevInMostPossible = stddev;
				mostPossibleStruct[0] = i;
				mostPossibleStruct[1] = j;
			}
		}
	}
	printf("Answer: %d, %d in %f.\n", mostPossibleStruct[0], mostPossibleStruct[1], stddevInMostPossible);
	
	return 0;
}

double printPeakToHKL(int* peaks, double* HKL, int peaksNum, int HKLNum){
	printf("%d, %d\n", peaksNum, HKLNum);
	for(int i = 0; i < peaksNum; i++){
		printf("%f, ", sin(peaks[i] * PI / 180.0));
	}
	printf("\n");
	if(HKLNum < peaksNum){
		peaksNum = HKLNum;
	}
	int HKLPeakIndex[peaksNum];	// Index of HKL correspond to each peaks.
	for(int i = 0; i < peaksNum; i++){
		HKLPeakIndex[i] = 0;
	}
	for(int i = 0; i < peaksNum - 1; i++){
		double peakDegreeRatio = sin(peaks[i + 1] * PI / 180.0) 
			/ sin(peaks[i] * PI / 180.0);
		// Possibility of HKL correspond to peak(lower then more possible).
		double minHKLIndex = fabs(1 - (HKL[HKLPeakIndex[i] + 1] / HKL[HKLPeakIndex[i]]) 
								  / peakDegreeRatio);
		HKLPeakIndex[i + 1] = HKLPeakIndex[i] + 1;
		for(int j = HKLPeakIndex[i] + 1; j < HKLNum; j++){
			double HKLIndex = fabs(1 - (HKL[j] / HKL[i]) / peakDegreeRatio);
			if(HKLIndex < minHKLIndex){
				minHKLIndex = HKLIndex;
				HKLPeakIndex[i + 1] = j;
			}
		}
	}
	
	double latticeConstAverage = 0;
	double latticeConstStddev = 0;
	for(int i = 0; i < peaksNum; i++){
		if(HKLPeakIndex[i] < peaksNum){
			printf("%d, (%d %f(in %d)), %f\n", i, peaks[i], HKL[HKLPeakIndex[i]], HKLPeakIndex[i
																							   ], sin(peaks[i] * PI / 180.0) / HKL[HKLPeakIndex[i]]);
			latticeConstAverage +=
				HKL[HKLPeakIndex[i]] * 0.154 / (sin(peaks[i] * PI / 180.0) * 2);
		}
	}
	latticeConstAverage /= peaksNum;
	for(int i = 0; i < peaksNum; i++){
		if(HKLPeakIndex[i] < peaksNum){
			double deferAverage = latticeConstAverage -
				HKL[HKLPeakIndex[i]] * 0.154 / (sin(peaks[i] * PI / 180.0) * 2);
			latticeConstStddev += deferAverage * deferAverage;
		}
	}
	latticeConstStddev /= peaksNum;
	printf("Average: %f, Stddev: %f\n", latticeConstAverage, sqrt(latticeConstStddev));
	return sqrt(latticeConstStddev);
}
