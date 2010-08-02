package com.newmacondo.georg.ICMM;

import java.util.HashSet;

public class CoreProfile {
	int mu;
	int delta[];
	int epsilon[][][];
	int gamma[][];
	ICMM myICMM = null;
	HashSet<Integer> genera = new HashSet<Integer>();
	
	public CoreProfile(ICMM mi) {
		myICMM = mi;
		delta = new int[myICMM.numTimepoints-1];
		epsilon = new int[myICMM.numGenera][myICMM.numTimepoints][myICMM.numReplicates];
		gamma = new int[myICMM.numGenera][myICMM.numTimepoints];
	}
}
