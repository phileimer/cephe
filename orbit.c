/* Jean Philippe EIMER - 1999 *
 *  phil.eimer@9online.fr     */

#include <math.h>

#include <stdio.h>

#include "include.h"


void	OrbiteInit(Orbite *o, Astre *pAst)
{
	static const double Vorb[9][6][4]={ \
		{{MER_A,0,0,0},MER_VE,MER_VI,MER_VL,MER_VO,MER_VOB}, \
		{{VEN_A,0,0,0},VEN_VE,VEN_VI,VEN_VL,VEN_VO,VEN_VOB}, \
		{{MAR_A,0,0,0},MAR_VE,MAR_VI,MAR_VL,MAR_VO,MAR_VOB}, \
		{{JUP_A,0,0,0},JUP_VE,JUP_VI,JUP_VL,JUP_VO,JUP_VOB}, \
		{{SAT_A,0,0,0},SAT_VE,SAT_VI,SAT_VL,SAT_VO,SAT_VOB}, \
		{{URA_A,0,0,0},URA_VE,URA_VI,URA_VL,URA_VO,URA_VOB}, \
		{{NEP_A,0,0,0},NEP_VE,NEP_VI,NEP_VL,NEP_VO,NEP_VOB}, \
		{{PLU_A,0,0,0},PLU_VE,PLU_VI,PLU_VL,PLU_VO,PLU_VOB}, \
		{{SOL_A,0,0,0},SOL_VE,SOL_VI,SOL_VL,SOL_VOB,SOL_VOB}};

	static const VectData	Vlun[5]={LUN_VO,LUN_VL,LUN_VM,LUN_VF,LUN_VD};
	static const VectData	Vlunper[4]={LUN_PER_VL,LUN_PER_VM,LUN_PER_VF,LUN_PER_VD};
	static const VectData	Vv[4]={LUN_PER_VA,LUN_PER_VB,LUN_PER_VC,LUN_PER_VE};

	
	o->pA = pAst;

	switch(o->pA->astreid)
	{
		case	LUN:
			// vecteurs de définition des éléments de l'orbite
			VectDataCopy(&o->VO, &Vlun[0]);
			VectDataCopy(&o->VL, &Vlun[1]);
			VectDataCopy(&o->VM, &Vlun[2]);
			VectDataCopy(&o->VF, &Vlun[3]);
			VectDataCopy(&o->VD, &Vlun[4]);

			// pour plus de précision
			VectDataCopy(&o->VLper, &Vlunper[0]);
			VectDataCopy(&o->VMper, &Vlunper[1]);
			VectDataCopy(&o->VFper, &Vlunper[2]);
			VectDataCopy(&o->VDper, &Vlunper[3]);
			VectDataCopy(&o->Va, &Vv[0]);
			VectDataCopy(&o->Vb, &Vv[1]);
			VectDataCopy(&o->Vc, &Vv[2]);
			VectDataCopy(&o->Ve, &Vv[3]);
			break;

		case	SOL:
			OrbiteCopyInit(o, Vorb[PLU+1]);
			break;

		case	TER:
			OrbiteCopy(o, &o->pA->pE->Soleil.ElemOrbit);
			o->L += 180.0;
			o->L = K360(o->L);
			o->Omega += 180.0;
			o->Omega = K360(o->Omega);
			o->omegab += 180.0;
			o->omegab = K360(o->omegab);
			break;

		default:
			OrbiteCopyInit(o, Vorb[o->pA->astreid]);
	}
}


// Initialisation des vecteurs
Orbite	*OrbiteCopyInit(Orbite *o, const double Vinit[6][4])
{
	o->a = Vinit[0][0];
	VectDataCopyInit(&o->Ve, Vinit[1]);
	VectDataCopyInit(&o->Vi, Vinit[2]);
	VectDataCopyInit(&o->VL, Vinit[3]);
	VectDataCopyInit(&o->VO, Vinit[4]);
	VectDataCopyInit(&o->Vob, Vinit[5]);

	return(o);
}


// Copie d'une structure Orbite
Orbite	*OrbiteCopy(Orbite *o, const Orbite *pO2)
{
	o->pA = pO2->pA;

	VectDataCopy(&o->VL, &pO2->VL);
	VectDataCopy(&o->VO, &pO2->VO);
	VectDataCopy(&o->Ve, &pO2->Ve);
	o->M = pO2->M;
	o->Omega = pO2->Omega;

	if (o->pA->astreid == LUN)
	{
		VectDataCopy(&o->VM, &pO2->VM);
		VectDataCopy(&o->VF, &pO2->VF);
		VectDataCopy(&o->VD, &pO2->VD);

		VectDataCopy(&o->VLper, &pO2->VLper);
		VectDataCopy(&o->VMper, &pO2->VMper);
		VectDataCopy(&o->VFper, &pO2->VFper);
		VectDataCopy(&o->VDper, &pO2->VDper);

		VectDataCopy(&o->Va, &pO2->Va);
		VectDataCopy(&o->Vb, &pO2->Vb);
		VectDataCopy(&o->Vc, &pO2->Vc);

		o->L = pO2->L;
		o->F = pO2->F;
		o->D = pO2->D;

		o->OmegaR = pO2->OmegaR;
		o->CCR = pO2->CCR;
		o->EE = pO2->EE;
	}
	else
	{
		VectDataCopy(&o->Vi, &pO2->Vi);
		VectDataCopy(&o->Vob, &pO2->Vob);

		o->a = pO2->a;
		o->e = pO2->e;
		o->i = pO2->i;
		o->L = pO2->L;
		o->omegab = pO2->omegab;
		o->deltaL = pO2->deltaL;

		o->aper = pO2->aper;
		o->eper = pO2->eper;
		o->Lper = pO2->Lper;
		o->vkper = pO2->vkper;
		o->rper = pO2->rper;
		o->lper = pO2->lper;
		o->bper = pO2->bper;

		o->MRsp = pO2->MRsp;
		o->omega = pO2->omega;

		o->E = pO2->E;
		o->ER = pO2->ER;

		o->r = pO2->r;
		o->v = pO2->v;
	}

	return(o);
}


// Calcul des éléments orbitaux en fonction du temps
void	OrbiteElements(Orbite *o, const Instant *pInst)
{
	if (o->pA->astreid == LUN)
	{
		// calcul des éléments de l'orbite
		o->Omega = VectDataDotTprecis2(&o->VO, pInst->t, pInst->datejulienne);
		o->OmegaR = DEG2RAD(o->Omega);
		o->L = VectDataDotTprecis2(&o->VL, pInst->t, pInst->datejulienne);
		o->L = K360(o->L);
		o->M = VectDataDotTprecis2(&o->VM, pInst->t, pInst->datejulienne);
		o->M = K360(o->M);
		o->F = VectDataDotTprecis2(&o->VF, pInst->t, pInst->datejulienne);
		o->F = K360(o->F);
		o->D = VectDataDotTprecis2(&o->VD, pInst->t, pInst->datejulienne);
		o->D = K360(o->D);

		o->CCR = o->OmegaR + DEG2RAD(VectDataDotT(&o->Vc, pInst->t));
		o->EE = VectDataDotT(&o->Ve, pInst->t);
	}
	else
	{
		// calcul à partir des vecteurs
		o->e = VectDataDotT(&o->Ve, pInst->t);
		o->i = VectDataDotT(&o->Vi, pInst->t);
		o->L = VectDataDotTprecis(&o->VL, pInst->t);
		o->Omega = VectDataDotT(&o->VO, pInst->t);
		o->omegab = VectDataDotT(&o->Vob, pInst->t);
		o->deltaL = o->VL.d[1]*360.0 + o->VL.d[2] + o->VL.d[3];

		o->L = K360(o->L);
		o->Omega = K360(o->Omega);
		o->omegab = K360(o->omegab);

		o->MRsp = DEG2RAD(o->L - o->omegab);	// anomalie moyenne sans perturbation, en radians
	}
}


void	OrbiteAnomalies(Orbite *o)
{
	// application des perturbations
	o->a += o->aper;
	o->omegab += o->vkper / o->e;
	o->e += o->eper;
	o->L += o->Lper;

	// calcul des autres éléments avec perturbations
	OrbiteAnomalieMoyenne(o);
	OrbiteKepler(o);
	OrbiteAnomalieVraie(o);
}


// calcul de M et omega
double	OrbiteAnomalieMoyenne(Orbite *o)
{
	o->M = o->L - o->omegab;
	o->omega = o->omegab - o->Omega;

	o->M = K360(o->M);
	o->omega = K360(o->omega);

	return(o->M);
}


// Equation de Kepler
double	OrbiteKepler(Orbite *o)
{
	double	MR=DEG2RAD(o->M);
	double	E2=o->ER;

	do
	{
		o->ER = E2;
		E2 = o->ER - (o->ER - o->e * sin(o->ER) - MR)/(1.0 - o->e * cos(o->ER));
	}
	while(fabs(E2 - o->ER) > PRECISION_KEPLER);

	o->E = RAD2DEG(E2);
	o->E = K360(o->E);

	o->ER = DEG2RAD(o->E);

	return(o->E);
}


// calcul du rayon vecteur et de l'anomalie vraie
double	OrbiteAnomalieVraie(Orbite *o)
{
	double	x,y;

	x = o->a * (cos(o->ER) - o->e);	// coordonnées rectangulaire
	y = o->a * sqrt(1 - o->e * o->e) * sin(o->ER);

	o->r = sqrt(x*x + y *y);	// coordonnées polaires
	o->v = RAD2DEG(atan(y/x));

	if(x<0.0)		// lever d'indétermination suite à atan
		o->v += 180.0; 	

	o->v = K360(o->v);		// 0-360°

	return(o->v);
}


// Calcul des perturbations sur les éléments orbitaux
void	OrbitePerturbations(Orbite *o, double t)
{
	o->aper = o->eper = o->Lper = o->vkper = o->rper = o->lper = o->bper = 0;	// Pas de perturbation par d�faut

	if(o->pA->astreid==SOL)
		OrbitePerturbationsSoleil(o, t);
	else if(o->pA->astreid < JUP)
		OrbitePerturbationsTellur(o, t);
	else if(o->pA->astreid < PLU)
		OrbitePerturbationsGeantes(o, t);
}

void	OrbitePerturbationsSoleil(Orbite *o, double t)
{
	static const VectData	Ve[4]=SOL_PER_VVE;
	static const VectData	Vper[2]={{SOL_PER_VL},{SOL_PER_VB}};
	static const double	Vbc[5]=SOL_PER_VBC;
	static double	Vlc[5]=SOL_PER_VLC;
	double	e[4];


	e[0]=DEG2RAD(VectDataDotTprecis(&Ve[0], t));
	e[1]=DEG2RAD(VectDataDotTprecis(&Ve[1], t));
	e[2]=DEG2RAD(VectDataDotTprecis(&Ve[2], t));
	e[3]=DEG2RAD(VectDataDotTprecis(&Ve[3], t));

	o->rper = Vbc[0]*sin(e[0])+Vbc[1]*sin(e[1])+Vbc[2]*sin(e[2])+Vbc[3]*cos(e[3])+Vbc[4]*sin(DEG2RAD(VectDataDotTprecis(&Vper[1], t)));
	o->lper = Vlc[0]*cos(e[0])+Vlc[1]*cos(e[1])+Vlc[2]*cos(e[2])+Vlc[3]*sin(e[3])+Vlc[4]*sin(DEG2RAD(VectDataDotT(&Vper[0], t)));
}

void	OrbitePerturbationsTellur(Orbite *o, double t)
{
	static const VectData	Vmer[8]=MER_PER_VV, Vven[12]=VEN_PER_VV, Vven2=VEN_PER_V;
	static const VectData	Vmar[21]=MAR_PER_VV, Vmar1=MAR_PER_V1, Vmar2=MAR_PER_V2;
	static const double	Cmer[8]=MER_PER_VC, Cven[12]=VEN_PER_VC;
	static const double	Cmar[21]=MAR_PER_VC1, Cmar2[2]=MAR_PER_VC2;
	double	aa;
	unsigned short int	i;
	VectData	VM;


	VM.d[0] = DEG2RAD(VectDataDotTprecis(&o->pA->pE->Planete[JUP].ElemOrbit.VL, t) - VectDataDotT(&o->pA->pE->Planete[JUP].ElemOrbit.Vob, t));	// anomalie moyenne de Jupiter
	VM.d[2] = o->MRsp;	// anomalie moyenne de Mercure si MER, de V�nus si VEN, de Mars si MAR
	VM.d[3] = 1.0;

	switch(o->pA->astreid)
	{
		case	MER:
			VM.d[1] = DEG2RAD(VectDataDotTprecis(&o->pA->pE->Planete[VEN].ElemOrbit.VL, t) - VectDataDotT(&o->pA->pE->Planete[VEN].ElemOrbit.Vob, t));	// anomalie moyenne de Venus
			for(i=0;i<4;i++)
				o->lper += Cmer[i] * cos(VectDataDot(&VM, &Vmer[i]));
			for(i=4;i<8;i++)
				o->rper += Cmer[i] * cos(VectDataDot(&VM, &Vmer[i]));
			break;

		case	VEN:
			VM.d[1] = DEG2RAD(VectDataDotTprecis(&o->pA->pE->Soleil.ElemOrbit.VL, t) - VectDataDotT(&o->pA->pE->Soleil.ElemOrbit.Vob, t));	// anomalie moyenne du  Soleil
			for(i=0;i<5;i++)
				o->lper += Cven[i] * cos(VectDataDot(&VM, &Vven[i]));
			for(i=5;i<12;i++)
				o->rper += Cven[i] * cos(VectDataDot(&VM, &Vven[i]));

			o->Lper = VEN_PER_C * sin(VectDataDotT(&Vven2, t));
			break;

		case	MAR:
			VM.d[1] = DEG2RAD(VectDataDotTprecis(&o->pA->pE->Soleil.ElemOrbit.VL, t) - VectDataDotT(&o->pA->pE->Soleil.ElemOrbit.Vob, t));	// anomalie moyenne du  Soleil
			for(i=0;i<8;i++)
				o->lper += Cmar[i] * cos(VectDataDot(&VM, &Vmar[i]));

			o->lper += MAR_PER_C1 * cos(VectDataDot(&VM, &Vmar1) - DEG2RAD(VectDataDotTprecis(&o->pA->pE->Planete[VEN].ElemOrbit.VL, t) - VectDataDotT(&o->pA->pE->Planete[VEN].ElemOrbit.Vob, t)));

			for(i=8;i<21;i++)
				o->rper += Cmar[i] * cos(VectDataDot(&VM, &Vmar[i]));

			aa = VectDataDot(&VM, &Vmar2);
			o->Lper = Cmar2[0] * sin(aa) + Cmar2[1] * cos(aa);
			break;
	}
}


void	OrbitePerturbationsGeantes(Orbite *o, double t)
{
	static const VectData	VJ[4]=GEA_PER_VJ;
	static const double	jL[33]=JUP_PER_VL, je[39]=JUP_PER_VE;
	static const double	jvk[22]=JUP_PER_VVK, ja[12]=JUP_PER_VA;
	static const double	sL[32]=SAT_PER_VL, se[68]=SAT_PER_VE;
	static const double	svk[27]=SAT_PER_VVK, sa[38]=SAT_PER_VA;
	static const double	sb[6]=SAT_PER_VB;
	static const double	uL[7]=URA_PER_VL, ue[5]=URA_PER_VE;
	static const double	uvk[4]=URA_PER_VVK, ur[13]=URA_PER_VR;
	static const double	ul[15]=URA_PER_VLH, ub[7]=URA_PER_VB;
	static const double	nL[5]=NEP_PER_VL, ne[5]=NEP_PER_VE, na[4]=NEP_PER_VA;
	static const double	nvk[4]=NEP_PER_VVK, nr[6]=NEP_PER_VR;
	static const double	nl[5]=NEP_PER_VLH, nb[2]=NEP_PER_VB;
	static double	j[14];
	static double	u[41];
	double	tmp;


	if(o->pA->astreid == JUP)	// premier passage = initialisation
	{
		j[0]=t/5.0+.1;		// j1
		j[1]=VectDataDotT(&VJ[0], t);	// j2
		j[2]=VectDataDotT(&VJ[1], t);
		j[3]=VectDataDotT(&VJ[2], t);	// j4

		j[4]=-2*j[1]+5*j[2];	// j5
		j[5]=2*j[1]-6*j[2]+3*j[3];	// j6

		j[6]=j[2]-j[1];		// j7

		j[7]=VectDataDotT(&VJ[3], t);	// j8 Uranus et Neptune
		j[8]=2*j[7]-j[3];	// j9 Ura Nept

		j[9]=j[3]-j[1];		// ja Uranus
		j[10]=j[3]-j[2];	// jb Uranus, j8 Saturne

		j[11]=j[7]-j[3];	// jc Uranus et Neptune
					
		j[12]=j[7]-j[1];	// ja Neptune
		j[13]=j[7]-j[2];	// jb Neptune

		//---------------------------------------------

		u[0]=sin(j[2]);		// u1
		u[1]=cos(j[2]);

		tmp=j[2]+j[2];
		u[2]=sin(tmp);		// u3
		u[3]=cos(tmp);

		u[4]=sin(j[4]);		// u5
		u[5]=cos(j[4]);
		u[6]=sin(j[4]+j[4]);	// u7
		u[7]=sin(j[5]);
		u[8]=sin(j[6]);		// u9
		u[9]=cos(j[6]);		// ua
		
		tmp=j[6]+j[6];
		u[10]=sin(tmp);		// ub
		u[11]=cos(tmp);		// uc

		tmp+=j[6];
		u[12]=sin(tmp);		// ud
		u[13]=cos(tmp);		// ue

		tmp+=j[6];
		u[14]=sin(tmp);		// uf
		u[15]=cos(tmp);		// ug

		u[16]=cos(tmp+j[6]);	// uh, vh

		tmp=3*j[2];
		u[17]=sin(tmp);		// ui
		u[18]=cos(tmp);		// uj

		tmp+=j[2];
		u[19]=sin(tmp);		// uk
		u[20]=cos(tmp);		// ul

		u[21]=cos(j[4]+j[4]);	// um, vi
		u[22]=sin(5*j[6]);	// un

		tmp=j[10]+j[10];
		u[23]=sin(tmp);		// uo
		u[24]=cos(tmp);		// up

		tmp+=j[10];
		u[25]=sin(tmp);		// uq
		u[26]=cos(tmp);		// ur

		u[27]=sin(j[8]);	// ut, vj
		u[28]=cos(j[8]);	// uu

		tmp=j[8]+j[8];
		u[29]=sin(tmp);		// uv
		u[30]=cos(tmp);		// uw

		u[31]=sin(j[10]);	// ux
		u[32]=cos(j[10]);	// uy
		u[33]=sin(j[3]);	// uz
		u[34]=cos(j[3]);	// va

		tmp=j[3]+j[3];
		u[35]=sin(tmp);		// vb
		u[36]=cos(tmp);		// vc

		tmp=j[11]+j[11];
		u[37]=sin(tmp);		// vd
		u[38]=cos(tmp);		// ve

		u[39]=sin(j[7]);	// vf
		u[40]=cos(j[7]);	// vg
	}

	switch(o->pA->astreid)
	{
		case	JUP:
			// perturbation en longitude moyenne (L)
			o->Lper = (jL[0]+(jL[1]+jL[2]*j[0])*j[0])*u[4];
			o->Lper += (jL[3]+(jL[4]+jL[5]*j[0])*j[0])*u[5];
			o->Lper += (jL[6]+(jL[7]+jL[8]*j[0])*j[0])*u[6];
			o->Lper += jL[9]*u[7]+jL[10]*u[8]+jL[11]*u[10];
			o->Lper += jL[12]*u[12]+jL[13]*u[14];
			tmp=jL[14]*u[10];
			tmp+=(jL[15]+jL[16]*j[0])*u[8]+jL[17]*u[12];
			tmp+=jL[20]*u[11];
			tmp+=(jL[21]+jL[22]*j[0])*u[9];
			o->Lper += tmp*u[0];
			tmp=(jL[18]+jL[19]*j[0])*u[8];
			tmp+=jL[23]*u[10];
			tmp+=(jL[24]*j[0]+jL[25])*u[9]+jL[26];
			tmp+=jL[27]*u[11]+jL[28]*u[13];
			o->Lper += tmp*u[1];
			o->Lper += (jL[29]*u[8]+jL[30]*u[10])*u[2];
			o->Lper += (jL[31]*u[9]+jL[32]*u[11])*u[3];

			// perturbation en excentricit� (e)
			o->eper = (je[0]+(je[1]+je[2]*j[0])*j[0])*u[4]+(je[3]+je[4]*j[0])*u[5];
			tmp=je[5]*u[8]+je[6]*u[10]+je[7]*u[12]+je[8];
			tmp+=(je[9]+je[10]*j[0])*u[9]+je[11]*u[11];
			o->eper += tmp*u[0];
			tmp=(je[12]+je[13]*j[0])*u[8]+je[14]*u[10]+je[15];
			tmp+=je[16]*u[9]+je[17]*u[11]+je[18]*u[13]+je[19]*u[15];
			tmp+=je[20]*u[16];
			o->eper += tmp*u[1];
			tmp=(je[21]+je[22]*j[0])*u[8]+je[23]*u[10];
			tmp+=je[24]*u[12]+(je[25]*j[0]+je[26])*u[9]+je[27]*u[11];
			tmp+=je[28]*u[13];
			o->eper += tmp*u[2];
			tmp=(je[29]*j[0]+je[30])*u[8]+je[31]*u[10];
			tmp+=je[32]*u[12]+je[33]+(je[34]+je[35]*j[0])*u[9];
			tmp+=je[36]*u[11]+je[37]*u[13];
			o->eper += tmp*u[3];
			o->eper *= je[38];

			// perturbation en longitude du périhélie (omegab)
			o->vkper = (jvk[0]+jvk[1]*j[0])*u[4];
			o->vkper += (j[0]*(jvk[3]*j[0]+jvk[4])+jvk[5])*u[5];
			tmp=jvk[2];
			tmp+=jvk[6]*u[9]+(jvk[7]+jvk[8]*j[0])*u[8];
			tmp+=jvk[9]*u[11]+jvk[10]*u[13];
			o->vkper += tmp*u[0];
			o->vkper += (jvk[11]*u[8]+jvk[12]*u[10]+jvk[13]*u[9])*u[1];
			o->vkper += (jvk[14]*u[8]+jvk[15]*u[10]+jvk[16]*u[9]+jvk[17]*u[11])*u[2];
			o->vkper += (jvk[18]*u[8]+jvk[19]*u[10]+jvk[20]*u[9]+jvk[21]*u[11])*u[3];

			// perturbation en demi grand axe (a)
			o->aper = ja[0]*u[9]+ja[1]*u[5]+ja[2]*u[11]+ja[3]*u[13]+ja[4]*u[15];
			o->aper += (ja[5]*u[8]+ja[6]*u[11])*u[0];
			o->aper += (ja[7]*u[10]+ja[8]*u[12]+ja[9]*u[9]+ja[10]*u[11])*u[1];
			o->aper *= ja[11];
			break;

		case	SAT:
			// perturbation en longitude moyenne (L)
			o->Lper = sL[0]*u[6]+sL[1]*u[7]+sL[2]*u[8];
			o->Lper += (sL[3]+(sL[4]+sL[5]*j[0])*j[0])*u[4];
			o->Lper += (sL[6]+(sL[7]+sL[8]*j[0])*j[0])*u[5];
			o->Lper += sL[9]*u[12]+sL[10]*u[14]+sL[11]*u[0];
			o->Lper += sL[13]*u[10];
			tmp=sL[12]*u[10];
			tmp+=(sL[14]+sL[15]*j[0])*u[8]+sL[16]*u[12];
			tmp+=(sL[17]+sL[18]*j[0])*u[9]+sL[19]*u[11];
			o->Lper += tmp*u[0];
			tmp=(sL[20]+sL[21]*j[0])*u[8]+sL[22]*u[11];
			tmp+=(sL[23]+sL[24]*j[0])*u[9]+sL[25]*u[13];
			o->Lper += tmp*u[1];
			o->Lper += (sL[26]*u[8]+sL[27]*u[10]+sL[28]*u[25])*u[2];
			o->Lper += (sL[29]*u[9]+sL[30]*u[11]+sL[31]*u[26])*u[3];

			// perturbation en excentricité (e)
			o->eper = (se[0]+(se[1]+se[2]*j[0])*j[0])*u[4];
			o->eper += (se[3]+(se[4]+se[5]*j[0])*j[0])*u[5]+(se[6]+se[7]*j[0])*u[6];
			o->eper += (se[8]+se[9]*j[0])*u[21]+se[10]*u[10];
			tmp=se[11]+(se[12]+se[13]*j[0])*u[8]+(se[14]+se[15]*j[0])*u[10];
			tmp+=se[16]*u[9]+se[17]*u[11]+se[18]*u[13]+se[19]*u[15];
			tmp+=se[20]*u[16]+se[21]*u[24];
			o->eper += tmp*u[0];
			tmp=se[22]+se[23]*j[0];
			tmp+=se[24]*u[8]+se[25]*u[10]+se[26]*u[12]+se[27]*u[14];
			tmp+=se[28]*u[22]+(se[29]+se[30]*j[0])*u[9];
			tmp+=(se[31]+se[32]*j[0])*u[11]+se[33]*u[23];
			o->eper += tmp*u[1];
			tmp=se[34]+(se[35]+se[36]*j[0])*u[8]+se[37]*u[10]+se[38]*u[12];
			tmp+=se[39]*u[14]+(se[40]+se[41]*j[0])*u[9];
			tmp+=(se[42]+se[43]*j[0])*u[11]+se[44]*u[13]+se[45]*u[25];
			tmp+=se[46]*u[26];
			o->eper += tmp*u[2];
			tmp=se[47]+(se[48]+se[49]*j[0])*u[8];
			tmp+=(se[50]+se[51]*j[0])*u[10]+se[52]*u[12];
			tmp+=(se[53]+se[54]*j[0])*u[9]+(se[55]+se[56]*j[0])*u[11];
			tmp+=se[57]*u[13]+se[58]*u[15]+se[59]*u[25]+se[60]*u[26];
			o->eper += tmp*u[3];
			o->eper += (se[61]*u[8]+se[62]*u[12])*u[17];
			o->eper += (se[63]*u[9]+se[64]*u[13])*u[18];
			o->eper += se[65]*u[13]*u[19]+se[66]*u[12]*u[20];
			o->eper *= se[67];
			
			// perturbation en longitude du périhélie (omegab)
			o->vkper = (svk[0]+(svk[1]+svk[2]*j[0])*j[0])*u[4];
			o->vkper += svk[3]*u[8];
			o->vkper += (svk[4]+(svk[5]+svk[6]*j[0])*j[0])*u[5];
			o->vkper += (svk[8]*u[8]+svk[9]*u[10]+svk[10]*u[12])*u[0];
			o->vkper += (svk[7]+svk[11]*u[9]+svk[12]*u[11]+svk[13]*u[13])*u[1];
			tmp=(svk[14]+svk[15]*j[0])*u[8];
			tmp+=(svk[17]+svk[18]*j[0])*u[9];
			tmp+=(svk[19]+svk[20]*j[0])*u[11];
			o->vkper += tmp*u[2];
			tmp=svk[16]*u[10]+(svk[21]+svk[22]*j[0])*u[8];
			tmp+=(svk[23]+svk[24]*j[0])*u[9];
			tmp+=(svk[25]+svk[26]*j[0])*u[11];
			o->vkper += tmp*u[3];

			// perturbation en demi grand axe (a)
			o->aper = sa[0]*u[4]+sa[2]*u[5]+sa[4]*u[9]+sa[6]*u[11];
			o->aper += sa[8]*u[13]+sa[11]*u[15]+sa[13]*u[16];
			tmp=sa[15]+sa[17]*u[8]+sa[19]*u[10]+sa[21]*u[12];
			tmp+=sa[23]*u[14]+sa[25]*u[9]+sa[27]*u[11];
			tmp+=sa[29]*u[13]+sa[31]*u[15];
			o->aper += tmp*u[0];
			tmp=sa[1]*u[10]+sa[3]*u[12]+sa[5]*u[14]+sa[7]*u[9];
			tmp+=(sa[9]+sa[10]*j[0])*u[11]+sa[12]*u[13];
			tmp+=sa[33]+sa[35]*u[8];
			o->aper += tmp*u[1];
			o->aper += (sa[14]*u[10]+sa[16]*u[9]+sa[18]*u[11]+sa[20]*u[13])*u[2];
			o->aper += (sa[22]*u[8]+sa[24]*u[10]+sa[26]*u[11]+sa[28]*u[13])*u[3];
			o->aper += (sa[30]*u[8]+sa[32]*u[12])*u[17];
			o->aper += (sa[34]*u[9]+sa[36]*u[13])*u[18];
			o->aper *= sa[37];

			// perturbation en latitude héliocentrique (b)
			o->bper = (sb[0]*u[0]+sb[1]*u[1])*u[9];
			o->bper += (sb[2]*u[10]+sb[3]*u[11])*u[2];
			o->bper += (sb[4]*u[10]+sb[5]*u[11])*u[3];
			break;

		case	URA:
			// perturbation en longitude moyenne (L)
			o->Lper = (uL[0]+uL[1]*j[0])*u[27];
			o->Lper += (uL[2]+uL[3]*j[0])*u[28]+uL[4]*u[29];
			o->Lper += uL[5]*u[30]+uL[6]*sin(j[5]);

			// perturbation en longitude du périhélie (omegab)
			o->vkper = uvk[0]*u[27]+uvk[1]*u[29];
			o->vkper += (uvk[2]+uvk[3]*j[0])*u[28];

			// perturbation en excentricité (e)
			o->eper = (ue[0]*j[0]+ue[1])*u[27]+ue[2]*u[28]+ue[3]*u[30];
			o->eper *= ue[4];

			// perturbation en demi grand axe (a)
			o->aper = URA_PER_A*u[28];

			// perturbation en longitude héliocentrique (l)
			o->lper = (ul[0]+(ul[1]+ul[2]*j[0])*j[0])*cos(j[3]+j[10]);
			o->lper += (ul[3]+ul[4]*j[0])*sin(j[3]+j[10]);
			o->lper += (ul[5]+(ul[6]+ul[7]*j[0])*j[0])*cos(2*j[3]+j[10]);
			o->lper += ul[8]*sin(j[3]+3*j[11])+ul[9]*sin(j[9]);
			o->lper += ul[10]*u[31]+ul[11]*u[32];
			o->lper += ul[12]*sin(j[11])+ul[13]*u[37];
			o->lper += ul[14]*sin(3*j[11]);

			// perturbation en latitude héliocentrique (b)
			o->bper = (ub[0]*u[31]+ub[1]*u[32]+ub[2]*cos(4*j[11]))*u[33];
			o->bper += (ub[3]*u[31]+ub[4]*u[32]+ub[5]*sin(4*j[10]))*u[34];
			o->bper += ub[6]*(u[38]*u[35]+u[37]*u[36]);
			
			// perturbation en rayon vecteur (r)
			o->rper = ur[0]+ur[1]*cos(j[9])+ur[2]*u[34]+ur[3]*u[32];
			o->rper += ur[4]*u[38]+ur[5]*(cos(j[11])-cos(3*j[11]));
			o->rper += (ur[6]*u[34]+ur[7]*u[33]+ur[8]*u[36])*u[31];
			o->rper += (ur[9]*u[34]+ur[10]*u[33]+ur[11]*u[35])*u[32];
			o->rper *= ur[12];
			break;

		case	NEP:
			// perturbation en longitude moyenne (L)
			o->Lper = (nL[0]*j[0]+nL[1])*u[27];
			o->Lper += (nL[2]*j[0]+nL[3])*u[28]+nL[4]*u[29];

			// perturbation en longitude du périhélie (omegab)
			o->vkper = nvk[0]*u[27]+nvk[1]*u[28]+nvk[2]*u[29];
			o->vkper += nvk[3]*u[30];

			// perturbation en excentricité (e)
			o->eper = ne[0]*u[27]+ne[1]*u[29]+ne[2]*u[28]+ne[3]*u[30];
			o->eper *= ne[4];

			// perturbation en demi grand axe (a)
			o->aper = na[0]*u[28]+na[1]*u[27]+na[2]*u[30];
			o->aper *= na[3];

			// perturbation en longitude héliocentrique (l)
			o->lper = nl[0]*sin(j[12])+nl[1]*sin(j[13]);
			o->lper += nl[2]*u[37]+nl[3]*u[38]*u[39]+nl[4]*u[37]*u[40];

			// perturbation en latitude héliocentrique (b)
			o->bper = nb[0]*u[38]*u[39]+nb[1]*u[37]*u[40];

			// perturbation en rayon vecteur (r)
			o->rper = nr[0]+nr[1]*cos(j[12])+nr[2]*cos(j[13]);
			o->rper += nr[3]*cos(j[11])+nr[4]*u[38];
			o->rper *= nr[5];
			break;
	}
}


// formatage des résultats
#if 0
void	OrbiteAffiche(Orbite *o)
{
	printf("\n%s", _("Eléments Orbitaux"));
	printf("\n%s\t\t a = %9.5f",_("Demi Grand Axe"), o->a);
	printf("\t%s\t\t e = %9.5f",_("Excentricité"), o->e);
	printf("\n%s\t i = %9.5f",_("Inclinaison / Eclipt."), o->i);

	printf("\t%s\t L = %9.5f",_("Longitude Moyenne"), o->L);
	printf("\n%s\t O = %9.5f",_("Longitude Noeud Ascend."), o->Omega);
	printf("\t%s\twb = %9.5f",_("Longitude Périhélie"), o->omegab);


	printf("\n%s\t M = %9.5f",_("Anomalie Moyenne"), o->M);
	printf("\t%s\t w = %9.5f",_("Argument Périhélie"), o->omega);

	printf("\n%s\t E = %9.5f",_("Anomalie Excentrique"), o->E);
	
	printf("\n%s\t\t r = %9.5f",_("Rayon Vecteur"), o->r);
	printf("\t%s\t\t v = %9.5f\n",_("Anomalie Vraie"), o->v);
}
#endif
