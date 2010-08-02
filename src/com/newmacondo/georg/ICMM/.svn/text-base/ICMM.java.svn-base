package com.newmacondo.georg.ICMM;

import java.io.IOException;
import java.util.ArrayList;

import umontreal.iro.lecuyer.probdist.GammaDist;

import com.newmacondo.georg.StatUtil.StatUtil;
import com.newmacondo.georg.Utils.FlatFileDataReader;

public class ICMM {
	int numGenera;
	int numTimepoints;
	int numReplicates;
	String[] generaNames;
	String[] timeLabels;
	ArrayList<CoreProfile> cores = new ArrayList<CoreProfile>();
	CoreProfile[] assignments = null;
	
	// we keep some cores around to avoid constant memory reallocation
	ArrayList<CoreProfile> bufferCores = new ArrayList<CoreProfile>();
	int minBufferCores = 10;
	int maxBufferCores = 25;
	
	// data Y_gtr (genera, timepoint, replicate)
	int[][][] Y;
	
	// prior mean for mu_k (empirical)
	double c_0 = 0.0;
	// prior variance for mu_k (empirical)
	double v_0 = 0.0;
	
	// variance for delta_{k t}
	double eta_delta[];
	// prior mean for eta_{delta t} (empirical variance)
	double v_delta[];
	// prior variance for eta_{delta t}
	double w_delta = 0.0;
	
	// variance for gamma_{g k t}
	double eta_gamma[];
	// prior mean for eta_{gamma t} (empirical variance)
	double v_gamma[];
	// prior variance for eta_{gamma t}
	double w_gamma = 0.0;
	
	// prior variance for replicates
	double v_R = 0.0;
	
	public void readData(String fileName) {
		FlatFileDataReader reader = new FlatFileDataReader();
		try {
			reader.readFileDiscrete(fileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		numTimepoints = reader.numCols;
		numReplicates = Integer.parseInt(reader.metaData.get("numReplicates"));
		numGenera = reader.numRows/numReplicates;
		Y = new int[numGenera][numTimepoints][numReplicates];
		generaNames = new String[numGenera];
		timeLabels = reader.colNames;
		
		int g = 0;
		int i = 0;
		int t = 0;
		int r = 0;
		
		for (i=0;i<reader.numRows;i++) {
			g = (int) Math.floor(((double) i)/((double) numReplicates));
			r = i % numReplicates;
			for (t=0;t<reader.numCols;t++) {
				Y[g][t][r] = reader.dvalues[i][t];
			}
			generaNames[g] = reader.rowNames[g];
		}
	}
	
	public void initHyperParameters() {
		int g = 0;
		int t = 0;
		int r = 0;
		int i = 0;
		double mn[][] = new double[numGenera][numTimepoints];
		double xt[] = new double[numGenera];
		
		eta_delta = new double[numTimepoints-1];
		v_delta = new double[numTimepoints-1];
		
		eta_gamma = new double[numTimepoints];
		v_gamma = new double[numTimepoints];
		
		for (g=0;g<numGenera;g++) {
			for (t=0;t<numTimepoints;t++) {
				for (r=0;r<numReplicates;r++) {
					mn[g][t] += (double) Y[g][t][r];
				}
				mn[g][t] = mn[g][t]/((double) numReplicates);
			}
		}
		
		// estimate c_0 and v_0
		for (g=0;g<numGenera;g++) {
			xt[g] = mn[g][0];
		}
		GammaDist tg = GammaDist.getInstanceFromMLE(xt,numGenera);
		c_0 = tg.getMean();
		v_0 = tg.getVariance();
		
		// estimate v_delta
		xt = new double[numGenera];
		for (t=1;t<numTimepoints;t++) {
			for (g=0;g<numGenera;g++) {
				if (mn[g][t-1] == 0) {
					xt[g] = mn[g][t];
				} else {
					xt[g] = mn[g][t]/mn[g][t-1];
				}
			}
			tg = GammaDist.getInstanceFromMLE(xt,numGenera);
			v_delta[t-1]  = tg.getVariance();
			eta_delta[t-1] = v_delta[t-1];
			System.out.println(v_delta[t-1]);
		}
		
		// set v_gamma = v_delta (v_gamma[0] is arbitrarily set to the median)
		v_gamma[0] = StatUtil.median(v_delta);
		eta_gamma[0] = v_gamma[0];
		for (t=1;t<numTimepoints;t++) {
			v_gamma[t] = v_delta[t-1];
			eta_gamma[t] = v_gamma[t];
		}
		
		// estimate v_R
		v_R = 0.0;
		for (g=0;g<numGenera;g++) {
			for (t=0;t<numTimepoints;t++) {
				for (r=0;r<numReplicates;r++) {
					if (Y[g][t][r] >  0) {
						v_R += Math.pow(1.0 - mn[g][t]/((double) Y[g][t][r]),2.0);
					} else {
						v_R += Math.pow(1.0-mn[g][t], 2.0);
					}
				}
			}
		}
		v_R = v_R /((double) numGenera*numTimepoints*numReplicates);
		System.out.println(v_R);
	}
}
