NEURON {
	POINT_PROCESS E2_Prob
	RANGE tau1, tau2, e, i, P, seed, rand, rel, total_rel
	NONSPECIFIC_CURRENT i
	RANGE g
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	e=0	(mV)
	seed=-99
	P=.24  :default value for P
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	rand
	rel
	total_rel
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
	if(seed>=0){	
		set_seed(seed)
	}
	rel = 0
	total_rel =0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
	rel=total_rel
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
	rand = scop_random()
	if (P>rand){	
		total_rel=total_rel+1			
		A = A + weight*factor
		B = B + weight*factor
	}
}
