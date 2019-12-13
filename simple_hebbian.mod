NEURON {
	POINT_PROCESS Hebb_E2Syn
	RANGE tau1, tau2, e, i, k, deltaw, interval, thresh
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
	thresh = 0 (mV)
	e=0	(mV)
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	tpost (ms)
	tpre (ms)
	interval (ms)
	deltaw
	factor
	k : learning rate
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
	tpost = -1e9
	deltaw = 0
	net_send(0, 1)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS), tpre (ms)) {
	INITIAL { tpre = -1e9 }
	if (flag == 0) { : presynaptic spike (after last post so depress)
		tpre = t
		interval = tpost-tpre
		deltaw = deltaw + k/interval
		A = A + (weight+deltaw)*factor
		B = B + (weight+deltaw)*factor	
	}else if (flag == 2) { : postsynaptic spike
		tpost = t
		FOR_NETCONS(w1, tp) {
			interval = t - tp
			deltaw = deltaw + k/interval
		}
	}else if (flag == 1) {
		WATCH (v > thresh) 2
	}
	

}
