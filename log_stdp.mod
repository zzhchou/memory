NEURON {
	POINT_PROCESS Log_E2Syn
	RANGE tau1, tau2, e, i, k, deltaw, interval, thresh
	RANGE alpha, b, c_p, c_d, tau_p, tau_d, f, J0
	RANGE seed
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
	alpha = 5
	b = 100
	c_p = 1
	c_n = 0.5
	tau_p = 17
	tau_n = 34
	J0 = 1e-6
	f
	seed = 0
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
	tpost = 0
	k = 0.1
	deltaw = 0
	set_seed(seed)
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

FUNCTION update(J, J0, interval) {
	if (interval > 0) {
		f = c_p*exp(-J/(J0*b))
		update = f*exp(-fabs(interval/tau_p))
	
	} else {
		if (J<=J0) {
			f = c_n*J/J0
		
		} else {
			f = c_n*(1+log(1+alpha*((J/J0)-1))/alpha)
		
		}
		
		update = -f*exp(-fabs(interval/tau_n))
	}
}

NET_RECEIVE(weight (uS), tpre (ms)) {
	INITIAL { tpre = -1e9 }
	if (flag == 0) { : presynaptic spike (after last post so depress)
		interval = tpost-t
		tpre = t
		if (tpost == 0) {
			deltaw = deltaw
		} else {
			deltaw = deltaw+k*(1+normrand(0,0.6))*update(deltaw+weight, J0, interval)
		}
		if (weight + deltaw <= 0) {
			deltaw = -weight
		}
		A = A + (weight+deltaw)*factor
		B = B + (weight+deltaw)*factor

	}else if (flag == 2) { : postsynaptic spike
		tpost = t
		FOR_NETCONS(w1, tp) {
			interval = t - tp
			tpost = t
			deltaw = deltaw+k*(1+normrand(0,0.6))*update(deltaw+w1, J0, interval)
		}
	}else if (flag == 1) {
		WATCH (v > thresh) 2
	}
	

}
