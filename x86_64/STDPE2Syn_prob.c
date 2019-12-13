/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__STDPE2_Prob
#define _nrn_initial _nrn_initial__STDPE2_Prob
#define nrn_cur _nrn_cur__STDPE2_Prob
#define _nrn_current _nrn_current__STDPE2_Prob
#define nrn_jacob _nrn_jacob__STDPE2_Prob
#define nrn_state _nrn_state__STDPE2_Prob
#define _net_receive _net_receive__STDPE2_Prob 
#define state state__STDPE2_Prob 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define tau1 _p[0]
#define tau2 _p[1]
#define e _p[2]
#define wmax _p[3]
#define wmin _p[4]
#define d _p[5]
#define p _p[6]
#define dtau _p[7]
#define ptau _p[8]
#define thresh _p[9]
#define gbdel _p[10]
#define gbint _p[11]
#define gblen _p[12]
#define d_count _p[13]
#define p_count _p[14]
#define A_init _p[15]
#define seed _p[16]
#define i _p[17]
#define g _p[18]
#define gs _p[19]
#define pot _p[20]
#define v_rec _p[21]
#define P _p[22]
#define rand _p[23]
#define rel _p[24]
#define total_rel _p[25]
#define C _p[26]
#define B _p[27]
#define tpost _p[28]
#define on _p[29]
#define factor _p[30]
#define DC _p[31]
#define DB _p[32]
#define _g _p[33]
#define _tsav _p[34]
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "d", 0, 1,
 "gbint", 1e-09, 1e+09,
 "gblen", 1e-09, 1e+09,
 "gbdel", 1e-09, 1e+09,
 "tau2", 1e-09, 1e+09,
 "tau1", 1e-09, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tau1", "ms",
 "tau2", "ms",
 "e", "mV",
 "wmax", "uS",
 "wmin", "uS",
 "dtau", "ms",
 "ptau", "ms",
 "thresh", "mV",
 "gbdel", "ms",
 "gbint", "ms",
 "gblen", "ms",
 "C", "uS",
 "B", "uS",
 "i", "nA",
 "g", "uS",
 0,0
};
 static double B0 = 0;
 static double C0 = 0;
 static double delta_t = 0.01;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
#define _watch_array _ppvar + 3 
 
#define _fnc_index 5
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   Prop* _prop = ((Point_process*)_vptr)->_prop;
   if (_prop) { _nrn_free_watch(_prop->dparam, 3, 2);}
   if (_prop) { _nrn_free_fornetcon(&(_prop->dparam[_fnc_index]._pvoid));}
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[6]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"STDPE2_Prob",
 "tau1",
 "tau2",
 "e",
 "wmax",
 "wmin",
 "d",
 "p",
 "dtau",
 "ptau",
 "thresh",
 "gbdel",
 "gbint",
 "gblen",
 "d_count",
 "p_count",
 "A_init",
 "seed",
 0,
 "i",
 "g",
 "gs",
 "pot",
 "v_rec",
 "P",
 "rand",
 "rel",
 "total_rel",
 0,
 "C",
 "B",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 35, _prop);
 	/*initialize range parameters*/
 	tau1 = 0.1;
 	tau2 = 10;
 	e = 0;
 	wmax = 0;
 	wmin = 0;
 	d = 0.5;
 	p = 0.5;
 	dtau = 34;
 	ptau = 17;
 	thresh = -20;
 	gbdel = 0;
 	gbint = 1000;
 	gblen = 1000;
 	d_count = 0;
 	p_count = 0;
 	A_init = 0;
 	seed = -99;
  }
 	_prop->param = _p;
 	_prop->param_size = 35;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 
#define _tqitem &(_ppvar[2]._pvoid)
 static void _net_receive(Point_process*, double*, double);
 extern int _nrn_netcon_args(void*, double***);
 static void _net_init(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _STDPE2Syn_prob_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 35, 7);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "watch");
  hoc_register_dparam_semantics(_mechtype, 4, "watch");
  hoc_register_dparam_semantics(_mechtype, 5, "fornetcon");
  hoc_register_dparam_semantics(_mechtype, 6, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 3;
 add_nrn_fornetcons(_mechtype, _fnc_index);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 STDPE2_Prob /home/zzchou/memory/x86_64/STDPE2Syn_prob.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   DC = - C / tau1 ;
   DB = - B / tau2 ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 DC = DC  / (1. - dt*( ( - 1.0 ) / tau1 )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau2 )) ;
  return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
    C = C + (1. - exp(dt*(( - 1.0 ) / tau1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - C) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tau2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - B) ;
   }
  return 0;
}
 
static double _watch1_cond(_pnt) Point_process* _pnt; {
  	_p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
	v = NODEV(_pnt->node);
	return  ( v ) - ( thresh ) ;
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{   int _watch_rm = 0;
    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   if ( _lflag  == 0.0 ) {
     rand = scop_random ( ) ;
     if ( P > rand ) {
         if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = C;
    double __primary = (C + ( _args[0] + _args[1] ) * factor) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau1 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - __primary );
    C += __primary;
  } else {
 C = C + ( _args[0] + _args[1] ) * factor ;
         }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B;
    double __primary = (B + ( _args[0] + _args[1] ) * factor) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau2 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - __primary );
    B += __primary;
  } else {
 B = B + ( _args[0] + _args[1] ) * factor ;
         }
 _args[2] = t ;
       if ( on  == 1.0 ) {
         _args[1] = _args[1] * ( 1.0 - d * exp ( ( tpost - t ) / dtau ) ) ;
         d_count = d_count + 1.0 ;
         }
       pot = _args[1] + _args[0] ;
       total_rel = total_rel + 1.0 ;
       }
     }
   else if ( _lflag  == 2.0  && on  == 1.0 ) {
     tpost = t ;
     {int _ifn1, _nfn1; double* _fnargs1, **_fnargslist1;
	_nfn1 = _nrn_netcon_args(_ppvar[_fnc_index]._pvoid, &_fnargslist1);
	for (_ifn1 = 0; _ifn1 < _nfn1; ++_ifn1) {
 	 _fnargs1 = _fnargslist1[_ifn1];
 {
       _fnargs1[1] = _fnargs1[1] + ( wmax - _fnargs1[0] - _fnargs1[1] ) * p * exp ( ( _fnargs1[2] - t ) / ptau ) ;
       pot = _fnargs1[1] + _fnargs1[0] ;
       p_count = p_count + 1.0 ;
       }
     	}}
 }
   else if ( _lflag  == 1.0 ) {
     {int _ifn2, _nfn2; double* _fnargs2, **_fnargslist2;
	_nfn2 = _nrn_netcon_args(_ppvar[_fnc_index]._pvoid, &_fnargslist2);
	for (_ifn2 = 0; _ifn2 < _nfn2; ++_ifn2) {
 	 _fnargs2 = _fnargslist2[_ifn2];
 {
       pot = _fnargs2[1] + _fnargs2[0] ;
       }
     	}}
   _nrn_watch_activate(_watch_array, _watch1_cond, 1, _pnt, _watch_rm++, 2.0);
 }
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
    _args[2] = - 1e9 ;
   }
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  B = B0;
  C = C0;
 {
   double _ltp ;
 if ( tau1 / tau2 > .9999 ) {
     tau1 = .9999 * tau2 ;
     }
   C = 0.0 ;
   B = 0.0 ;
   _ltp = ( tau1 * tau2 ) / ( tau2 - tau1 ) * log ( tau2 / tau1 ) ;
   factor = - exp ( - _ltp / tau1 ) + exp ( - _ltp / tau2 ) ;
   factor = 1.0 / factor ;
   gs = 1.0 ;
   on = 1.0 ;
   tpost = - 1e9 ;
   pot = 0.0 ;
   net_send ( _tqitem, (double*)0, _ppvar[1]._pvoid, t +  0.0 , 1.0 ) ;
   v_rec = v ;
   P = .24 ;
   if ( seed >= 0.0 ) {
     set_seed ( seed ) ;
     }
   rel = 0.0 ;
   total_rel = 0.0 ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = B - C ;
   i = g * gs * ( v - e ) ;
   gs = gs ;
   d_count = d_count ;
   p_count = p_count ;
   pot = pot ;
   v_rec = v ;
   rel = total_rel ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 { error =  state();
 if(error){fprintf(stderr,"at line 90 in file STDPE2Syn_prob.mod:\n	SOLVE state METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 }}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(C) - _p;  _dlist1[0] = &(DC) - _p;
 _slist1[1] = &(B) - _p;  _dlist1[1] = &(DB) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/zzchou/memory/STDPE2Syn_prob.mod";
static const char* nmodl_file_text = 
  ": STDP by Hines, changed to dual exponential (BPG 6-1-09)\n"
  ": Modified by BPG 13-12-08\n"
  ": Limited weights: max weight is wmax and min weight is wmin\n"
  ": (initial weight is specified by netconn - usually set to wmin)\n"
  ": Rhythmic GABAB suppresses conductance and promotes plasticity.\n"
  ": When GABAB is low, conductance is high and plasticity is off.\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS STDPE2_Prob\n"
  "	RANGE tau1, tau2, e, i, d, p, dtau, ptau, thresh, wmax, wmin\n"
  "	RANGE g, gbdel, gblen, gbint, B, C, gs, v_rec, d_count, p_count, pot, A_init,rel, total_rel\n"
  "	RANGE P, rand, seed\n"
  "	NONSPECIFIC_CURRENT i\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(uS) = (microsiemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	tau1=.1 (ms) <1e-9,1e9>\n"
  "	tau2 = 10 (ms) <1e-9,1e9>\n"
  "	e = 0	(mV)\n"
  "	wmax = 0 (uS)\n"
  "	wmin = 0 (uS)	: not used - use netconn weight instead (BPG)\n"
  "	d = 0.5 <0,1>: depression factor (multiplicative to prevent < 0)\n"
  "	p = 0.5 : potentiation factor (additive, non-saturating)\n"
  "	dtau = 34 (ms) : depression effectiveness time constant\n"
  "	ptau = 17 (ms) : Bi & Poo (1998, 2001)\n"
  "	thresh = -20 (mV)	: postsynaptic voltage threshold\n"
  "	gbdel = 0 (ms) <1e-9,1e9> : initial GABAB off interval (ms)\n"
  "	gbint = 1000 (ms) <1e-9,1e9> : GABAB off interval (ms)\n"
  "	gblen = 1000 (ms) <1e-9,1e9> : GABAB on length (ms)\n"
  "	d_count=0\n"
  "	p_count=0\n"
  "	A_init=0\n"
  "	seed=-99\n"
  "	\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	i (nA)\n"
  "	tpost (ms)\n"
  "	on\n"
  "	g (uS)\n"
  "	gs\n"
  "	factor\n"
  "	pot\n"
  "	v_rec\n"
  "	P\n"
  "	rand\n"
  "	rel\n"
  "	total_rel\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	C (uS)\n"
  "	B (uS)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	LOCAL tp\n"
  "	if (tau1/tau2 > .9999) {\n"
  "		tau1 = .9999*tau2\n"
  "	}\n"
  "	C = 0\n"
  "	B = 0\n"
  "	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)\n"
  "	factor = -exp(-tp/tau1) + exp(-tp/tau2)\n"
  "	factor = 1/factor\n"
  "	gs=1\n"
  "	on=1	: initially plastic\n"
  "	tpost = -1e9\n"
  "	pot=0\n"
  "	net_send(0, 1)\n"
  "	v_rec=v\n"
  "	P=.24\n"
  "	if(seed>=0){	\n"
  "		set_seed(seed)\n"
  "	}\n"
  "	rel = 0\n"
  "	total_rel =0\n"
  "		\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	g = B - C\n"
  "	i = g*gs*(v - e)\n"
  "	gs=gs\n"
  "	d_count=d_count\n"
  "	p_count=p_count\n"
  "	pot=pot\n"
  "	v_rec=v\n"
  "	rel=total_rel\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "	C' = -C/tau1\n"
  "	B' = -B/tau2\n"
  "}\n"
  "\n"
  "NET_RECEIVE(w (uS), A, tpre (ms)) {\n"
  "	INITIAL { tpre = -1e9 }\n"
  "	if (flag == 0) { : presynaptic spike  (after last post so depress)\n"
  "		rand = scop_random()\n"
  "		if (P > rand) {\n"
  ":			printf(\"entry flag=%g t=%g w=%g A=%g tpre=%g tpost=%g\\n\", flag, t, w, A, tpre, tpost)\n"
  "	:		g = g + w + A	: only for single exp (BPG)\n"
  "			C = C + (w + A)*factor\n"
  "			B = B + (w + A)*factor\n"
  "			tpre = t\n"
  "			if (on == 1) {\n"
  "				A = A * (1 - d*exp((tpost - t)/dtau))\n"
  "				d_count=d_count+1\n"
  "			}\n"
  "			pot=A+w\n"
  "			total_rel=total_rel+1\n"
  "		}\n"
  "	}else if (flag == 2 && on == 1) { : postsynaptic spike\n"
  ":		printf(\"entry flag=%g t=%g tpost=%g\\n\", flag, t, tpost)\n"
  "		tpost = t\n"
  "		FOR_NETCONS(w1, A1, tp) { : also can hide NET_RECEIVE args\n"
  ":			printf(\"entry FOR_NETCONS w1=%g A1=%g tp=%g  \", w1, A1, tp)\n"
  "			A1 = A1 + (wmax-w1-A1)*p*exp((tp - t)/ptau)\n"
  ":			printf(\"A1 increased to %g \\n\", A1)\n"
  "			pot=A1+w1\n"
  "			p_count=p_count+1\n"
  "		}\n"
  "	} else if (flag == 1) { : flag == 1 from INITIAL block\n"
  ":		printf(\"entry flag=%g t=%g\\n\", flag, t)\n"
  "		FOR_NETCONS(w1, A1, tp) { : also can hide NET_RECEIVE args\n"
  "			pot=A1+w1\n"
  "		}\n"
  "		WATCH (v > thresh) 2\n"
  "	} \n"
  "	\n"
  "}\n"
  ;
#endif
