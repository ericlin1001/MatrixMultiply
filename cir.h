//* There are three parts in the program
//* 1)cir.cpp is the circuit simulation in time domain using Billy method.
//*		input:RC,C2,RC2,C3,R2;	output: Vout[31252];
//*	2)test.cpp is calculated values of overshoot and settlingtime.
//*		inout: Vout[31252];		output: overshoot, settlingtime;
//* 3)cir_ga.cpp is main program which include GA program.
//*		input: overshoot, settlingtime;		output:RC,C2,RC2,C3,R2, fitness value;
//*		ZHANG JUN , 1999, City University of Hong Kong
//* Using Dual-loop Fuzzy Tuning GA seach method to find the best optimal parameter set.
//* The finial Version
//* ZHANG JUN , 1999 , City University of Hong Kong
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <conio.h>
#include <time.h>
#include <string.h>
#include <iostream>
using namespace std;
//*******************************begin of cir.h********************************
void cir();
int  chk_stage();
void do_stage();
void endofcycle();
void EnforceSWc();
double Eref();
void feedback();
void findend();
double findramp(int,double);
double findSWinstant(int, double);
void imgUpdateSW();
void LCbranch(double *);
void output();
void pop_value();
void pqa();
void push_value();
void restart_stage();
void setinit();
void simulateTo(double);
void snapOneSWcycle();
void snapThrough(double);
void Spqa();
void stage1();
void stage2();
void stage3();
void stage4();
void UpdateK();
void UpdateSW();
double Vg();
//*******************************end of cir.hpp********************************

//*****************************begin of cir()*************************************
#define h1 0e-3
#define h2 30e-3
#define h3 0.8*h2
#define zero 1e-13
#define nX 4
#define nV 6
#define nSW 2
#define delta 1e-6
#define mf 3

//#define Vg1 20	
//#define R 5.0

double Vg1;
double R;

#define RE 0.005
#define RD 0.05	
#define RT 0.05	

//#define L 0.00020000 //0.000138000//200e-6
//#define C 0.001000000 //0.000513000//1e-3

#define RL 	0.25
#define RC	0.1 
//#define RC3	x_ga[1] //4.7e3
//#define C2	x_ga[2]	//2e-6
#define RC2	/*x_ga[3]	//*/0.1
//#define C3	x_ga[3]//4]	//3.3e-6
//#define R2	x_ga[4]//5]	//300e3

//#define C4 x_ga[5]//6]	//1.8e-6
//#define R4 x_ga[6]//7]	//1000


#define RC4 0.1

//#define R1 x_ga[7]//0.6e3
double ACCEPT = 150;
#define Eref1 5.0
#define alpha 0.2e6
#define Imax 4.0
#define Ts 50e-6
#define Dt 0.85*Ts
FILE *output1;                    /* X[0] output */
FILE *output2;                    /* Vout output*/
int stage, Vcon, Denforce, condition,/*=9,*/ k;
double times, t;
double time_sim[32000];
double RC3,C2,/*RC2,*/C3,R2;
double C4,R4,R1;
double L=0.00020000, C=0.001000000;
double Vout[31252];
double X[nX], V[nV], SW[nSW];
double X0[nX],V0[nV],SW0[nSW];
double X2[nX],V2[nV],SW2[nSW];
double K0[nX], K1[nX];

double	overshoot,undershoot,PF1,OF4,OF2,ov,uv,time_con;
int settlingtime=31252;
double settling_time=h2;
double p_function;
///////////////////////////////////////////

double STD_R[]={100,110,120,130,150,160,180,200,220,240,270,300,330,360,390,430,470,510,560,620,680,750,820,910,
	1E3,1.1E3,1.2E3,1.3E3,1.5E3,1.6E3,1.8E3,2E3,2.2E3,2.4E3,2.7E3,3E3,3.2E3,3.3E3,3.6E3,3.9E3,4.3E3,4.7E3,5.1E3,5.6E3,6.2E3,6.6E3,7.5E3,8.2E3,9.1E3,
	10E3,11E3,12E3,13E3,15E3,16E3,18E3,20E3,22E3,24E3,27E3,30E3,33E3,36E3,39E3,43E3,47E3,51E3,56E3,62E3,68E3,75E3,82E3,91E3,100E3};
/*double STD_R[]={62,68,75,82,91,100,110,120,130,150,160,180,200,220,240,270,300,330,360,390,430,470,510,560,620,680,750,820,910,
  1E3,1.1E3,1.2E3,1.3E3,1.5E3,1.6E3,1.8E3,2E3,2.2E3,2.4E3,2.7E3,3E3,3.2E3,3.3E3,3.6E3,3.9E3,4.3E3,4.7E3,5.1E3,5.6E3,6.2E3,6.6E3,7.5E3,8.2E3,9.1E3,
  10E3,11E3,12E3,13E3,15E3,16E3,18E3,20E3,22E3,24E3,27E3,30E3,33E3,36E3,39E3,43E3,47E3,51E3,56E3,62E3,68E3,75E3,82E3,91E3,100E3,
  110E3,120E3,130E3,150E3,160E3,180E3,200E3,220E3,240E3,270E3,300E3,330E3,360E3,390E3,430E3,470E3,510E3,560E3,620E3,680E3,750E3,820E3,910E3,
  1E6,1.1E6,1.2E6,1.3E6,1.5E6,1.6E6,1.8E6,2E6,2.2E6,2.4E6,2.7E6,3E6};*/
double STD_C[]={1.0,1.1,1.2,1.3,1.5,1.6,1.8,2.0,2.2,2.4,2.7,3.0,3.3,3.6,3.9,4.3,4.7,5.1,5.6,6.2,6.8,7.5,8.2,9.1,10};
//STD_C[]*10^n
int RS=sizeof(STD_R)/sizeof(double);//RS=74
int CS=sizeof(STD_C)/sizeof(double);//CS=24;
//*****************************end of cir()******************************
//=======================test()===================================
//Test is used to calculate the overshoot, settling time, undershoot and so on.
void test()
{
	int	i,i_max;  // set the value to be zero at the first time
	int	sse_count=0;	// count the steady error in order to look for oscillation
	int	num=31252;
	int	str_count=0;	// start transient coun

	//double	v[31252],vv[31252],sse,s_under_s,v_filter[31252];
	double	v[32000],vv[32000],v_filter[32000];
	double E2=0.0;
	double area_ripple=0.0;
	double	s_over_s=0;		// start over shoot
	double	up_limit,lo_limit,max,min;
	double v_ref=5.0; //5.0783;//which is calculaed by c:\aa\ave.m
	double err[32000],err_sum;
	memset(err,0,32000*sizeof(double));//BUG FIXED::err should be initialized to 0.

	p_function=0;

	err_sum=0;

	up_limit=v_ref*(1+0.02); // reference band
	lo_limit=v_ref*(1-0.02);// reference band

	i_max=0;
	max=0.0;//up_limit;
	min=0.0;//lo_limit;	// max and min value to detect over & under shoot

	//settling_time=h2;
	//settlingtime=num;
	undershoot=0.0;
	overshoot=0.0;

	for (i=2;i<num;i=i+1) {v[i]=Vout[i];}   //Let v[i] = Vout

	for (i=5000;i<30000;i=i+1)	//Calculate error area
	{
		err[i]=fabs(v[i]-v_ref);
		err_sum=err_sum+err[i];
	}

	for (i=0;i<30000;i=i+1)
	{
		if((v[i]>up_limit)||(v[i]<lo_limit))
		{ 
			str_count=i;	//find settling time
			//		settling_time=time_sim[i];
			settling_time=str_count*(1e-6);
		}
	}

	for (i=0;i<30000;i=i+1)
	{
		if (v[i]>(max+v_ref)) 
		{
			max=v[i]-v_ref;  //find overshoot			
			i_max=i;	//find overshoot point
		}

		s_over_s=max;
	}

	for (i=0;i<30000;i=i+1)
	{
		if (((v_ref-v[i])>=min)&&(i>i_max)&&(i_max>1))
		{
			min=v_ref-v[i];
		}
	}

	overshoot=s_over_s;
	settlingtime=str_count;//point of settling time
	settling_time=settling_time;//real time , unit is sec
	undershoot=min;

	//Objective Function 1  OF1 is min steady state error  filter problem
	//The steady state error of v[i] within the required input voltage range 
	//and output load range
	//v_filrer[i] is the output after filter	
	v_filter[0]= v[4999];
	vv[0]=v[4999];
	vv[1]=v[5000];
	v_filter[1]= 0.98*v_filter[0]+ 0.5*(1-0.98)*(vv[1] + vv[0]); 
	for (i = 2;i<(num-5000);i=i+1)
	{
		vv[i]=v[(i+5000)];
		v_filter[i] = 0.98*v_filter[i-1]+ 0.5*(1-0.98)*(vv[i]+ vv[i-1]); 
		E2=E2+(v_filter[i]-v_ref)*(v_filter[i]-v_ref);
	}
	PF1=2.0*exp(-E2/500.);

	//Objective Function 2  OF2 is min overshoot, undershoot and settling time
	// The maximum overhoot and undershoot , damping ration and the settling time
	// of v[i]during the start up of the circuit
	ov=4.0/(1+exp((0.1-overshoot/v_ref)/(-0.013647)));
	uv=4.0/(1+exp((0.05-undershoot/v_ref)/(-0.006833)));
	OF2=4.0/(1+exp(-(0.005-settling_time)/0.0008));	//OF2 is for settling time

	//Objective Function 3  OF3 is min ripple
	for (i=5000;i<num;i=i+1)
	{
		if (err[i]>(0.02*v_ref))
			area_ripple=area_ripple+err[i]-(0.02*v_ref);
	}
	OF4=2.0*exp(-area_ripple/500.0);
}
//============end of test()=============================================================

//=============cir()====================================================================
//double Vout[31252];
//   V[0]=supply voltage, V[1]=diode voltage, V[2]=output voltage,
//   V[3]=op-amp inputs, V[4]=op-amp output, V[5]=ramp function
//   X[0]=iL, X[1]=Vc, X[2]=Vc2,X[3]=Vc3,SW[0]=iT, SW[1]=iD.
void output()
{
	/*if ((output1 = fopen("out1.dat","a"))==NULL)
	  {exit(1);}
	  if ((output2 = fopen("out2.dat","a"))==NULL)
	  {exit(1);}
	  if (times>h1) 
	  {fprintf(output1,"%f  %f\n",times,X[0]);
	  fprintf(output2, "%f\n",V[2]);  */
	if (k<31252)
	{
		time_sim[k]=times;
		Vout[k]=V[2];
		k=k+1;
	}
	/*  }
		fclose(output1);
		fclose(output2); */
}
double Vg()
{
	double VG=Vg1;
	return VG;
}
double Eref()
{
	double ER= Eref1*(1-exp(-times/((R4+RC4)*C4)));
	return ER;
}
void push_value()
{
	int i;
	for (i=0;i<nV;i++)   V0[i]=V[i];
	for (i=0;i<nSW;i++) SW0[i]=SW[i];
	for (i=0;i<nX;i++)   X0[i]=X[i];
}
void pop_value()
{
	int i;
	for (i=0;i<nV;i++)   V[i]=V0[i];
	for (i=0;i<nSW;i++) SW[i]=SW0[i];
	for (i=0;i<nX;i++)   X[i]=X0[i];
}
void pqa()
{
	for (int i=0;i<nX;i++)
		X[i]+=delta*(3.*K0[i]-K1[i])/2.;
}
void Spqa()
{
	for (int i=0;i<nX;i++)
		X[i]+=delta*(3.*K0[i]-K1[i])/2./mf;
}
void LCbranch(double *ptr)
{  //definition of L,C branch
	ptr[0]=(V[1]-V[2]-RL*X[0])/L;
	ptr[1]=(V[2]-X[1])/RC/C;
	ptr[2]=(V[3]-V[4]-X[2])/RC2/C2;
	ptr[3]=(V[3]-V[4]-X[3])/RC3/C3;  
}
void UpdateK()
{
	for (int i=0;i<nX;i++) K1[i]=K0[i];
	LCbranch(K0);
}
void feedback()
{
	V[3]=Eref();
	V[2]=(X[0]+X[1]/RC+V[3]/R1)/(1/RC+1/R+1/R1);
	V[4]=(V[3]*(1./R2+1./R1+1./RC2+1./RC3)-V[2]/R1-X[2]/RC2-X[3]/RC3)/(1./R2+1./RC2+1./RC3);
	V[5]=alpha*t;
}
void stage1()
{ // T on , D off
	V[0]=Vg()-RE*X[0];
	V[1]=V[0]-RT*X[0];
	feedback();
}
void stage2()
{ // T off, D on
	V[0]=Vg();
	V[1]=-RD*X[0];
	feedback();
}
void stage4()
{ // T on, D on
	V[1]=(Vg()-(RT+RE)*X[0])/(RT+RE-RD);
	V[0]=V[1]+RT*(V[1]/RD+X[0]);
	feedback();
}
void stage3()
{ // T off, D off
	V[0]=Vg();
	feedback();
	V[1]=RL*X[0]+V[2];
}
void do_stage()
{
	switch (stage)
	{
		case 1: stage1(); break;
		case 2: stage2(); break;
		case 3: stage3(); break;
		case 4: stage4(); break;
		default: cout<<"Undefined stage happened in do_stage!\n";
				 cin.get(); exit(1);
	}
}
void UpdateSW()
{
	int Vcon0=Vcon;
	Vcon= ((!Denforce)&&(t<Dt)&&(V[5]<V[4]))? 1:0;
	if (Vcon)
	{
		SW[0]=(V[0]-V[1])/RT;
		SW[1]=-V[1]/RD;
	}
	else
	{
		SW[0]=(V[0]-V[1])/RT;
		SW[1]=X[0];
		if (SW[1]==0.) SW[1]=-V[1]/RD;
	}
	if ((Vcon0)&&(SW[0]>Imax)) {Denforce=1; Vcon=0;}
}
void imgUpdateSW()
{
	SW[0]=(V[0]-V[1])/RT; //MOSFET
	SW[1]=-V[1]/RD; //diode
}
int chk_stage()
{
	int S,Ton,Don;
	Ton= (SW[0]>=0.)? 1: 0;
	Don= (SW[1]>0.)?  1: 0;
	if ((Vcon)&&(Ton)) S= (Don)? 4:1;
	else S= (Don)? 2:3;
	return S;
}
void EnforceSWc()
{ //switching condition enforcement
	switch(condition)
	{
		case 1: X[0]=SW[1]=0.; break;
		case 2: SW[0]=-zero; break;
		case 3: t=Dt; break;
		case 4: V[5]=V[4]+zero; UpdateSW(); break;
		case 5: Denforce=1; break;
		case 9: t=0.; Denforce=0; break;
		default: cout<<"No condition fulfills in EnforceSWc!\n";
				 cin.get(); exit(1);
	}
}
void restart_stage()
{
	double sV0[nX],sV1[nX];
	int newstage, i;
	stage=chk_stage();
	do_stage();
	UpdateSW();
	EnforceSWc();
	newstage=chk_stage();
	while (stage!=newstage)
	{
		stage=newstage;
		do_stage(); UpdateSW();
		EnforceSWc();
		newstage=chk_stage();
	}
	// the following is the algorithm fo find d/dt(x(-1))
	do_stage(); UpdateSW(); EnforceSWc();
	push_value();
	LCbranch(sV0); //first point.
	for (i=0;i<nX;i++) K1[i]=K0[i]=sV0[i]; //step fitting
	Spqa(); do_stage(); LCbranch(K0);
	for (i=1;i<mf;i++)
	{ Spqa(); do_stage(); UpdateK();}
	LCbranch(sV1);
	for (i=0;i<nX;i++)
	{
		K0[i]=sV0[i];
		K1[i]=2.*sV0[i]-sV1[i];
	}
	pop_value();
}
void snapThrough(double STEP)
{ // considering from K0, X0, X2
	double kh1,kh2;
	kh1=STEP*delta*(1.-STEP); kh2=STEP*STEP;
	for (int i=0;i<nX;i++) X[i]=K0[i]*kh1+kh2*(X2[i]-X0[i])+X0[i];
}
double findSWinstant(int SWp, double ref)
{
	double a,b,c,sol,B,D;
	a=(SW0[SWp]-2.*SW[SWp]+SW2[SWp])*2.;
	b=(-3.*SW0[SWp]+4.*SW[SWp]-SW2[SWp]);
	c=SW0[SWp]-ref;
	if (a!=0.)
	{ B=-b/2./a;
		D=B*B-c/a;
		if (D>=0.) D=sqrt(D);
		else  { cout<<"Negative determinant exists in findSW!\n";
			cin.get(); exit(1); }
		sol= (B>D)? B-D:B+D;
	}
	else //linear cord
		sol= -c/b;
	return sol;
}
double findramp(int Vp, double gnd)
{
	double a,b,c,B,D,sol;
	a=(V0[Vp]-2.*V[Vp]+V2[Vp])*2.;
	b=(-3.*V0[Vp]+4.*V[Vp]-V2[Vp])-alpha*delta;
	c=V0[Vp]-gnd;
	if (a!=0.)
	{ B=-b/2./a;
		D=B*B-c/a;
		if (D>=0.) D=sqrt(D);
		else  { cout<<"error happened in findramp!\n";
			cin.get(); exit(1); }
		sol= (B>D)? B-D: B+D;
	}
	else //linear cord
		sol= -c/b;
	return sol;
}
void findend()
{
	double STEP=2., newSTEP;
	int i;
	t-=delta; times-=delta;          //return back
	for (i=0;i<nX;i++) X2[i]=X[i];  //push up t(n+1) to t(n+2)
	for (i=0;i<nV;i++) V2[i]=V[i];
	for (i=0;i<nSW;i++) SW2[i]=SW[i];
	//finding values in the middle of such step.
	for (i=0;i<nX;i++) X[i]=(X2[i]+K0[i]*delta+X0[i]*3.)*0.25;
	times+=delta/2.; t+=delta/2.;
	do_stage(); imgUpdateSW();
	times-=delta/2.; t-=delta/2.;
	if ((SW0[1]>0.)&&(SW2[1]<0.)) //diode current becomes zero.
	{ STEP=findSWinstant(1,0.); condition=1;}
	if ((SW0[0]>0.)&&(SW2[0]<0.)) //transistor current becomes zero.
	{ newSTEP=findSWinstant(0,0.);
		if (newSTEP<STEP) {STEP=newSTEP; condition=2;} }
	if ((t<Dt)&&(t>=Dt-delta))  //maximum duty cycle limitation
	{ newSTEP=(Dt-t)/delta;
		if (newSTEP<STEP) {STEP=newSTEP; condition=3;} }
	if ((V0[4]>=V0[5])&&(V2[4]<V2[5])) //+ve ramp cutting
	{ newSTEP=findramp(4,V0[5]);
		if (newSTEP<STEP) {STEP=newSTEP; condition=4;} }
	if ((SW0[0]<Imax)&&(SW2[0]>Imax)) //maximum current
	{ newSTEP=findSWinstant(0,Imax);
		if (newSTEP<STEP) {STEP=newSTEP; condition=5;} }
	if (!condition)  //none of the above conditions happen!
	{ cout<<"No valid boundary condition fulfill in findend!\n";
		cin.get(); exit(1); }
	snapThrough(STEP);
	times+=STEP*delta; t+=STEP*delta;
	do_stage(); UpdateSW();
	restart_stage();
	output();
}
void endofcycle()
{
	double LEN=Ts-t;
	times+=LEN;
	for (int i=0;i<nX;i++) X2[i]=X[i];
	snapThrough(1.+LEN/delta);
}
void setinit()
{
	for (int i=0;i<nX;i++) X[i]=0.;
}
void snapOneSWcycle()
{
	int WithinCycle=1;
	condition=9; EnforceSWc();
	restart_stage();
	while (WithinCycle)
	{
		times+=delta; t+=delta;
		push_value();
		pqa(); do_stage(); UpdateSW();
		if (stage!=chk_stage()) findend();
		else if (t>Ts) {endofcycle(); WithinCycle=0;}
		else UpdateK();
		output();
	}
}
void simulateTo(double endtime)
{
	setinit(); output();
	while ((times<endtime)) snapOneSWcycle();
}
void cir()
{
	k=0;
	times=0;
	condition=9;
	simulateTo(h2);
}
//============end of cir()==============================================================
//============================================================================================
void init_parameters()
{    
	//================================
	int i;
	Vg1 = R = RC3 = C2 = C3 = R2 = C4 = R4 = R1 = 0.;
	stage = Vcon = Denforce = condition = k = 0;
	times = t = 0;
	for(i = 0; i < 32000; i++) time_sim[i] = 0.0;
	for(i = 0; i < 31252; i++) Vout[i] = 0.0;
	for(i = 0; i < nX; i++)
		X[i] = V[i] = X0[i] = V0[i] = 
			X2[i] = V2[i] = K0[i] = K1[i] = 0.0;
	for(i = 0; i < nSW; i++)
		SW[i] = SW0[i] = SW2[i] = 0.0;

	overshoot = undershoot = PF1 = OF4 = 
		OF2 = ov = uv = time_con = 0;
	settlingtime=31252;
	settling_time=h2;
	p_function = 0;
}

double _objective(double x[],double y[])
{	
	double fitness = 0.0;
	init_parameters();
	double u;
	int i,j;
	for(i=0;i<4;i++){//resistor
		for(j=0;j<RS;j++){
			if(x[i]<STD_R[j])
				break;
		}
		if(j==0) y[i]=STD_R[j];
		else if(j==RS) y[i]=STD_R[j-1];
		else if((x[i]-STD_R[j-1])<(STD_R[j]-x[i])) y[i]=STD_R[j-1];
		else y[i]=STD_R[j];
	}
	int cnt;
	for(i=4;i<7;i++){//capacitor
		cnt=0;
		u=x[i];
		while(u<1.0){u*=10;cnt++;}
		for(j=0;j<CS;j++){
			if(u<STD_C[j])
				break;
		}
		if(j==0) y[i]=STD_C[j];
		else if(j==CS) y[i]=STD_C[j-1];
		else if((u-STD_C[j-1])<(STD_C[j]-u)) y[i]=STD_C[j-1];
		else y[i]=STD_C[j];
		y[i]=y[i]*pow(10,-cnt);
	}
	R1= y[0];
	R2= y[1];
	RC3= y[2];
	R4= y[3];
	C2= y[4];
	C3= y[5];
	C4= y[6];
	///////////////
	time_con=2.0*exp(-(R4*C4/ (( (1.8e-6)*1000.0)*5)));	  	
	for (R=2;R<12;R=R+3)
	{ 
		for(Vg1=20;Vg1<61;Vg1=Vg1+20)
		{

			cir();
			test();

			/*			 if(settlingtime<20000)
						 {
						 time_con=2.0*exp(-(R4*C4/ (( (1.8e-6)*1000.0)*5)));	  
						 }
						 else
						 {
						 time_con=0.0;
						 }
			 */
			fitness +=PF1+uv+ov+OF2+OF4+time_con;
			//			 fprintf(out_ga1, "\n%5d,      %6.9f,%6.9f, %6.9f, %9d, %6.9f ", generation,  R, Vg1, overshoot,settlingtime,settling_time);
			//			 fprintf(out_ga1, "\n%6.9f, %6.9f, %6.9f, %6.9f, %6.9f ",PF1,uv,ov,OF2,OF4);
		}
	}
	return fitness;
}
double _objective(double x[])
{	
	double fitness = 0.0;
	init_parameters();
	R1= x[0];
	R2= x[1];
	RC3= x[2];
	R4= x[3];
	C2= x[4];
	C3= x[5];
	C4= x[6];
	///////////////
	time_con=2.0*exp(-(R4*C4/ (( (1.8e-6)*1000.0)*5)));	  	
	for (R=2;R<12;R=R+3)
	{ 
		for(Vg1=20;Vg1<61;Vg1=Vg1+20)
		{

			cir();
			test();

			/*			 if(settlingtime<20000)
						 {
						 time_con=2.0*exp(-(R4*C4/ (( (1.8e-6)*1000.0)*5)));	  
						 }
						 else
						 {
						 time_con=0.0;
						 }
			 */
			fitness +=PF1+uv+ov+OF2+OF4+time_con;
			//			 fprintf(out_ga1, "\n%5d,      %6.9f,%6.9f, %6.9f, %9d, %6.9f ", generation,  R, Vg1, overshoot,settlingtime,settling_time);
			//			 fprintf(out_ga1, "\n%6.9f, %6.9f, %6.9f, %6.9f, %6.9f ",PF1,uv,ov,OF2,OF4);
		}
	}
	return fitness;	
}
