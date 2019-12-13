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
 
#define nrn_init _nrn_init__Mult_E2Syn
#define _nrn_initial _nrn_initial__Mult_E2Syn
#define nrn_cur _nrn_cur__Mult_E2Syn
#define _nrn_current _nrn_current__Mult_E2Syn
#define nrn_jacob _nrn_jacob__Mult_E2Syn
#define nrn_state _nrn_state__Mult_E2Syn
#define _net_receive _net_receive__Mult_E2Syn 
#define state state__Mult_E2Syn 
 
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
#define thresh _p[2]
#define e _p[3]
#define alpha _p[4]
#define b _p[5]
#define c_p _p[6]
#define tau_p _p[7]
#define J0 _p[8]
#define f _p[9]
#define seed _p[10]
#define i _p[11]
#define g _p[12]
#define interval _p[13]
#define deltaw _p[14]
#define k _p[15]
#define A _p[16]
#define B _p[17]
#define tpost _p[18]
#define tpre _p[19]
#define factor _p[20]
#define DA _p[21]
#define DB _p[22]
#define _g _p[23]
#define _tsav _p[24]
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
 static double _hoc_update();
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
 "update", _hoc_update,
 0, 0
};
#define update update_Mult_E2Syn
 extern double update( double , double , double );
 /* declare global and static user variables */
#define c_n c_n_Mult_E2Syn
 double c_n = 1.65;
#define tau_n tau_n_Mult_E2Syn
 double tau_n = 34;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "tau2", 1e-09, 1e+09,
 "tau1", 1e-09, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tau1", "ms",
 "tau2", "ms",
 "thresh", "mV",
 "e", "mV",
 "A", "uS",
 "B", "uS",
 "i", "nA",
 "g", "uS",
 "interval", "ms",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double delta_t = 0.01;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "c_n_Mult_E2Syn", &c_n_Mult_E2Syn,
 "tau_n_Mult_E2Syn", &tau_n_Mult_E2Syn,
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
"Mult_E2Syn",
 "tau1",
 "tau2",
 "thresh",
 "e",
 "alpha",
 "b",
 "c_p",
 "tau_p",
 "J0",
 "f",
 "seed",
 0,
 "i",
 "g",
 "interval",
 "deltaw",
 "k",
 0,
 "A",
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
 	_p = nrn_prop_data_alloc(_mechtype, 25, _prop);
 	/*initialize range parameters*/
 	tau1 = 0.1;
 	tau2 = 10;
 	thresh = 0;
 	e = 0;
 	alpha = 5;
 	b = 100;
 	c_p = 1;
 	tau_p = 17;
 	J0 = 0.25;
 	f = 0;
 	seed = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 25;
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

 void _mult_stdp_reg() {
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
  hoc_register_prop_size(_mechtype, 25, 7);
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
 pnt_receive_size[_mechtype] = 2;
 add_nrn_fornetcons(_mechtype, _fnc_index);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Mult_E2Syn /home/zzchou/memory/x86_64/mult_stdp.mod\n");
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
   DA = - A / tau1 ;
   DB = - B / tau2 ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 DA = DA  / (1. - dt*( ( - 1.0 ) / tau1 )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau2 )) ;
  return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
    A = A + (1. - exp(dt*(( - 1.0 ) / tau1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - A) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tau2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - B) ;
   }
  return 0;
}
 
double update (  double _lJ , double _lJ0 , double _linterval ) {
   double _lupdate;
 if ( _linterval > 0.0 ) {
     f = c_p ;
     _lupdate = f * exp ( - fabs ( _linterval / tau_p ) ) ;
     }
   else {
     f = c_n * _lJ ;
     _lupdate = - f * exp ( - fabs ( _linterval / tau_n ) ) ;
     }
   
return _lupdate;
 }
 
static double _hoc_update(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  update (  *getarg(1) , *getarg(2) , *getarg(3) );
 return(_r);
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
     interval = tpost - t ;
     _args[1] = t ;
     if ( tpost  == 0.0 ) {
       deltaw = deltaw ;
       }
     else {
       deltaw = deltaw + k * ( 1.0 + normrand ( 0.0 , 0.6 ) ) * update ( _threadargscomma_ deltaw + _args[0] , J0 , interval ) ;
       }
     if ( _args[0] + deltaw <= 0.0 ) {
       deltaw = - _args[0] ;
       }
       if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A;
    double __primary = (A + ( _args[0] + deltaw ) * factor) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau1 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - __primary );
    A += __primary;
  } else {
 A = A + ( _args[0] + deltaw ) * factor ;
       }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B;
    double __primary = (B + ( _args[0] + deltaw ) * factor) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau2 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - __primary );
    B += __primary;
  } else {
 B = B + ( _args[0] + deltaw ) * factor ;
       }
 }
   else if ( _lflag  == 2.0 ) {
     tpost = t ;
     {int _ifn1, _nfn1; double* _fnargs1, **_fnargslist1;
	_nfn1 = _nrn_netcon_args(_ppvar[_fnc_index]._pvoid, &_fnargslist1);
	for (_ifn1 = 0; _ifn1 < _nfn1; ++_ifn1) {
 	 _fnargs1 = _fnargslist1[_ifn1];
 {
       interval = t - _fnargs1[1] ;
       tpost = t ;
       deltaw = deltaw + k * ( 1.0 + normrand ( 0.0 , 0.6 ) ) * update ( _threadargscomma_ deltaw + _fnargs1[0] , J0 , interval ) ;
       }
     	}}
 }
   else if ( _lflag  == 1.0 ) {
       _nrn_watch_activate(_watch_array, _watch1_cond, 1, _pnt, _watch_rm++, 2.0);
 }
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
    _args[1] = - 1e9 ;
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
  A = A0;
  B = B0;
 {
   double _ltp ;
 if ( tau1 / tau2 > .9999 ) {
     tau1 = .9999 * tau2 ;
     }
   A = 0.0 ;
   B = 0.0 ;
   _ltp = ( tau1 * tau2 ) / ( tau2 - tau1 ) * log ( tau2 / tau1 ) ;
   factor = - exp ( - _ltp / tau1 ) + exp ( - _ltp / tau2 ) ;
   factor = 1.0 / factor ;
   tpost = 0.0 ;
   k = 0.1 ;
   deltaw = 0.0 ;
   set_seed ( seed ) ;
   net_send ( _tqitem, (double*)0, _ppvar[1]._pvoid, t +  0.0 , 1.0 ) ;
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
   g = B - A ;
   i = g * ( v - e ) ;
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
 if(error){fprintf(stderr,"at line 67 in file mult_stdp.mod:\n	SOLVE state METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 }}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(A) - _p;  _dlist1[0] = &(DA) - _p;
 _slist1[1] = &(B) - _p;  _dlist1[1] = &(DB) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/zzchou/memory/mult_stdp.mod";
static const char* nmodl_file_text = 
  "NEURON {\n"
  "	POINT_PROCESS Mult_E2Syn\n"
  "	RANGE tau1, tau2, e, i, k, deltaw, interval, thresh\n"
  "	RANGE alpha, b, c_p, c_d, tau_p, tau_d, f, J0\n"
  "	RANGE seed\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	RANGE g\n"
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
  "	thresh = 0 (mV)\n"
  "	e=0	(mV)\n"
  "	alpha = 5\n"
  "	b = 100\n"
  "	c_p = 1\n"
  "	c_n = 1.65\n"
  "	tau_p = 17\n"
  "	tau_n = 34\n"
  "	J0 = 0.25\n"
  "	f\n"
  "	seed = 0\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	i (nA)\n"
  "	g (uS)\n"
  "	tpost (ms)\n"
  "	tpre (ms)\n"
  "	interval (ms)\n"
  "	deltaw\n"
  "	factor\n"
  "	k : learning rate\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	A (uS)\n"
  "	B (uS)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	LOCAL tp\n"
  "	if (tau1/tau2 > .9999) {\n"
  "		tau1 = .9999*tau2\n"
  "	}\n"
  "	A = 0\n"
  "	B = 0\n"
  "	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)\n"
  "	factor = -exp(-tp/tau1) + exp(-tp/tau2)\n"
  "	factor = 1/factor\n"
  "	tpost = 0\n"
  "	k = 0.1\n"
  "	deltaw = 0\n"
  "	set_seed(seed)\n"
  "	net_send(0, 1)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	g = B - A\n"
  "	i = g*(v - e)\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "	A' = -A/tau1\n"
  "	B' = -B/tau2\n"
  "}\n"
  "\n"
  "FUNCTION update(J, J0, interval) {\n"
  "	if (interval > 0) {\n"
  "		f = c_p\n"
  "		update = f*exp(-fabs(interval/tau_p))	\n"
  "	} else {\n"
  "		f = c_n*J\n"
  "		update = -f*exp(-fabs(interval/tau_n))\n"
  "	}\n"
  "}\n"
  "\n"
  "NET_RECEIVE(weight (uS), tpre (ms)) {\n"
  "	INITIAL { tpre = -1e9 }\n"
  "	if (flag == 0) { : presynaptic spike (after last post so depress)\n"
  "		interval = tpost-t\n"
  "		tpre = t\n"
  "		if (tpost == 0) {\n"
  "			deltaw = deltaw\n"
  "		} else {\n"
  "			deltaw = deltaw+k*(1+normrand(0,0.6))*update(deltaw+weight, J0, interval)\n"
  "		}\n"
  "		if (weight + deltaw <= 0) {\n"
  "			deltaw = -weight\n"
  "		}\n"
  "		A = A + (weight+deltaw)*factor\n"
  "		B = B + (weight+deltaw)*factor\n"
  "\n"
  "	}else if (flag == 2) { : postsynaptic spike\n"
  "		tpost = t\n"
  "		FOR_NETCONS(w1, tp) {\n"
  "			interval = t - tp\n"
  "			tpost = t\n"
  "			deltaw = deltaw+k*(1+normrand(0,0.6))*update(deltaw+w1, J0, interval)\n"
  "		}\n"
  "	}else if (flag == 1) {\n"
  "		WATCH (v > thresh) 2\n"
  "	}\n"
  "	\n"
  "\n"
  "}\n"
  ;
#endif
