/*
	A simple code for itebd, which aims to recover some results of Vidal's original work,
	namely, fig5&6 of arXiv:cond-mat/0605597

	This code is based on the ITensor Library v3.
*/

#include "itensor/all.h"
#include "fstream"
#include "iostream"
#include "math.h"

using namespace itensor;
using namespace std;

ITensor inverM (ITensor diaM);
ITensor copyTensor(ITensor term, IndexSet indnew); //only valid for diagTensor
ITensor GSAB(ITensor ga, ITensor la, ITensor gb, ITensor lb);
ITensor GSBA(ITensor gb, ITensor lb, ITensor ga, ITensor la);
ITensor gateU(int roi, int aob, Real field, Real tstep, SiteSet sitesh);
Real measureOp(ITensor ga, ITensor la, ITensor gb, ITensor lb, string oper, SiteSet sites);
void upDateA(ITensor& ga, ITensor& la, ITensor& gb, const ITensor& lb, ITensor gate, int bondset);
void upDateB(ITensor& gb, ITensor& lb, ITensor& ga, const ITensor& la, ITensor gate, int bondset);


// NOTATIONS:
//#########################################################################################
/*
	(lb)--LB--(la')--GA(sa)--(la)--LA--(lb')--GB(sb)--(lb)--LB--(la') for GSAB
	(la)--LA--(lb')--GB(sb)--(lb)--LB--(la')--GA(sa)--(la)--LA--(lb') for GSBA
	sa & sb is generated on SiteSet HalfSpin(2)
*/
//#########################################################################################


int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		printfln("Usage: %s input_file", argv[0]);
		return 0;
	}

	cout << "The program is running!!!\n";

	auto input = InputGroup(argv[1], "input");
	int d = 2;
	auto h0 = input.getReal("h0"); //field for fig5
	auto h1 = input.getReal("h1"); //initial field for fig6
	auto h2 = input.getReal("h2"); //changed field for fig6
	int precise1 = input.getInt("precise1"); //GS control for fig5
	int loop1 = input.getInt("loop1");
	int precise2 = input.getInt("precise2"); //GS control for gif6
	int loop2 = input.getInt("loop2");
	int bond = input.getInt("bond"); //bond for links
	auto tstep = input.getReal("tstep");
	auto ttotal = input.getReal("ttotal");
	int length = input.getInt("length"); //how long the two-point correlation is calculated

	//define the inds used
	auto sites = SpinHalf(2, {"ConserveQNs=",false});
	auto sa = sites(1);
	auto sb = sites(2);
	auto la = Index(bond, "LinkA");
	auto lb = Index(bond, "LinkB");

/************************** initialization of GA, GB, LA, LB **************************/
	auto LB = delta({lb, prime(la)});
	LB /= sqrt(bond);
	auto AB = randomITensor({sa, lb, sb, prime(la)});
	auto [U, S, V] = svd(AB, {sa,lb}, {"MaxDim",bond});
	auto lbond = commonIndex(U, S);
	auto rbond = commonIndex(S, V);
	auto LA = copyTensor(S, IndexSet(la, prime(lb)));
	LA /= norm(LA);
	auto inLB = inverM(LB);
	U = U*delta(lbond, la);
	V = V*delta(rbond, prime(lb));
	auto GA = U*inLB;
	auto GB = V*inLB;

/**************************************  FIG 6  ****************************************/
	ofstream outfile;
	outfile.open("fig6.dat");

	ITensor gateA, gateB;
	for(int digit = 0; digit <= precise2; ++digit)
	{
		gateA = gateU(2, 1, h1, pow(10, -digit), sites);
		gateB = gateU(2, 2, h1, pow(10, -digit), sites);
		for(int iter = 1; iter <= loop2; ++iter)
		{
			upDateA(GA, LA, GB, LB, gateA, bond);
			upDateB(GB, LB, GA, LA, gateB, bond);
		}
		cout << "loop of h1 GS: " << digit << "\n";
	}

	auto psiA = LB*GA*LA;
	auto psiB = LA*GB*LB;
	Real mz = measureOp(GA, LA, GB, LB, "Sz", sites);
	outfile << 0.0 << "\t" << abs(mz) << "\n";
	cout << "h1 ground state is found! Time simulation starts!!!" << "\n";

	//external field is changed to be h2
	gateA = gateU(1, 1, h2, tstep/2, sites);
	gateB = gateU(1, 2, h2, tstep, sites);
	for(int iter = 1; iter <= int(ttotal/tstep); ++iter)
	{
		cout << "time step is: " << tstep*iter << "\n";
		upDateA(GA, LA, GB, LB, gateA, bond);
		upDateB(GB, LB, GA, LA, gateB, bond);
		upDateA(GA, LA, GB, LB, gateA, bond);
		mz = measureOp(GA, LA, GB, LB, "Sz", sites);
		outfile << tstep*iter << "\t" << abs(mz) << "\n";
	}
	cout << "fig6 is done!!!" << "\n";
	outfile.close();

/**************************************  FIG 5  ****************************************/
	outfile.open("fig5.dat");
	cout << "fig5's calculation starts!!!" << "\n";

	for(int digit = 0; digit <= precise1; ++digit)
	{
		gateA = gateU(2, 1, h0, pow(10, -digit), sites);
		gateB = gateU(2, 2, h0, pow(10, -digit), sites);
		for(int iter = 1; iter <= loop1; ++iter)
		{
			upDateA(GA, LA, GB, LB, gateA, bond);
			upDateB(GB, LB, GA, LA, gateB, bond);
		}
		cout << "loop of h0 GS: " << digit << "\n";
	}
	cout << "GS is found!!!" << "\n";

	// two point correlation function is calculated iteratively
	auto m0 = measureOp(GA, LA, GB, LB, "Sz", sites);

	ITensor proxy, swap, RAZ, RBZ;
	proxy = GA*LA;
	RAZ = proxy*op(sites,"Sz",1);
	RAZ = noPrime(RAZ,"Site")*dag(noPrime(proxy,"LinkA"));
	proxy = GB*LB;
	RBZ = proxy*op(sites,"Sz",2);
	RBZ = noPrime(RBZ,"Site")*dag(noPrime(proxy,"LinkB"));

	proxy = LB*GA*LA;
	swap = noPrime(proxy*op(sites,"Sz",1),"Site");
	proxy = swap*dag(prime(proxy,"LinkB,1"));
	proxy.noPrime("LinkB,2");
	for(int iter = 1; iter <= length; ++iter)
	{
		if(iter%2 == 1)
		{
			mz = realPart(proxy*RBZ).elt();
			swap = GB*LB;
			proxy = (proxy*swap)*dag(noPrime(noPrime(swap,"LinkA"),"LinkB"));
		}
		else
		{
			mz = realPart(proxy*RAZ).elt();
			swap = GA*LA;
			proxy = (proxy*swap)*dag(noPrime(noPrime(swap,"LinkA"),"LinkB"));
		}
		outfile << iter << "\t" << 4.*(mz-pow(m0,2)) << "\n";
		cout << iter << "th correlation;\n";
	}
	cout << "finished!!!" << "\n";
	outfile.close();

	return 0;
}


//#######################################################################################
/*------------------------------define the functions used------------------------------*/
//#######################################################################################


ITensor GSAB(ITensor ga, ITensor la, ITensor gb, ITensor lb)
{
	// function to compute the ground state of configuration LB--GA--LA--GB--LB
	auto proxy1 = la*gb*lb;
	auto proxy2 = lb*ga;
	auto gs = proxy1*proxy2;
	return gs;
}

ITensor GSBA(ITensor gb, ITensor lb, ITensor ga, ITensor la)
{
	// function to compute the ground state of configuration LA--GB--LB--GA--LA
	auto proxy1 = lb*ga*la;
	auto proxy2 = la*gb;
	auto gs = proxy1*proxy2;
	return gs;
}

Real measureOp(ITensor ga, ITensor la, ITensor gb, ITensor lb, string oper, SiteSet sites)
{
	// function to compute the expectation value of "oper"
	ITensor proxy, swap;
	proxy = lb*ga*la;
	swap = proxy*op(sites, oper, 1);
	swap *= dag(prime(proxy,"Site"));
	auto m = realPart(swap).elt();
	proxy = la*gb*lb;
	swap = proxy*op(sites, oper, 2);
	swap *= dag(prime(proxy,"Site"));
	m += realPart(swap).elt();
	return m/2.;
}

ITensor copyTensor(ITensor term, IndexSet indnew)
{
	// diagTensor can not be contracted with delta, which is a sad story :(
	auto indold = inds(term);
	auto newterm = ITensor(indnew);
	for(int comp = 1; comp <= maxDim(term); ++comp)
	{
		newterm.set(indnew(1)=comp, indnew(2)=comp, term.elt(indold(1)=comp, indold(2)=comp));
	}
	return newterm;
}

ITensor inverM(ITensor diaM)
{
	// function to get the inverse of LA and LB
	auto inv = [](Real x) { if(x == 0) return 0.; else return 1./x; };
	diaM.apply(inv);
	return diaM;
}

ITensor gateU(int roi, int aob, Real field, Real tstep, SiteSet sitesh)
{
	// function to generate the gate U for RealTime or ImagTime; for AB or BA configuration
	// if the model is TI by two sites, AB & BA should be treated differently as below
	ITensor result, hterm;
	if(aob == 1)
	{
		hterm = 4*op(sitesh,"Sx",1)*op(sitesh,"Sx",2);
		hterm += field*( op(sitesh,"Sz",1)*op(sitesh,"Id",2) + op(sitesh,"Id",1)*op(sitesh,"Sz",2) );
		if(roi == 1) result = BondGate(sitesh, 1, 2, BondGate::tReal, tstep, hterm);
		else result = BondGate(sitesh, 1, 2, BondGate::tImag, tstep, hterm);
	}
	else
	{
		hterm = 4*op(sitesh,"Sx",2)*op(sitesh,"Sx",1);
		hterm += field*( op(sitesh,"Sz",2)*op(sitesh,"Id",1) + op(sitesh,"Id",2)*op(sitesh,"Sz",1) );
		if(roi == 1) result = BondGate(sitesh, 2, 1, BondGate::tReal, tstep, hterm);
		else result = BondGate(sitesh, 2, 1, BondGate::tImag, tstep, hterm);
	}
	return result;
}

void upDateA(ITensor& ga, ITensor& la, ITensor& gb, const ITensor& lb, ITensor gate, int bondset)
{
	// function to update for LB--GA--LA--GB--LB configuration
	auto gsab = GSAB(ga, la, gb, lb);
	gsab *= gate;
	gsab.noPrime("Site");
	auto sa = findIndex(ga,"Site");
	auto ba = commonIndex(ga, la);
	auto bb = commonIndex(gb, lb);
	auto [U, S, V] = svd(gsab, {sa,bb}, {"MaxDim",bondset});
	auto lbond = commonIndex(U, S);
	auto rbond = commonIndex(S, V);
	auto ldelta = delta(lbond, ba);
	auto rdelta = delta(rbond, prime(bb));
	S = copyTensor(S, IndexSet(ba, prime(bb)));
	U *= ldelta;
	V *= rdelta;
	la = S/norm(S);
	auto inlb = inverM(lb);
	ga = U*inlb;
	gb = V*inlb;
}

void upDateB(ITensor& gb, ITensor& lb, ITensor& ga, const ITensor& la, ITensor gate, int bondset)
{
	// function to update for LA--GB--LB--GA--LA configuration
	auto gsba = GSBA(gb, lb, ga, la);
	gsba *= gate;
	gsba.noPrime("Site");
	auto sb = findIndex(gb,"Site");
	auto ba = commonIndex(ga, la);
	auto bb = commonIndex(gb, lb);
	auto [U, S, V] = svd(gsba, {sb, ba}, {"MaxDim",bondset});
	auto lbond = commonIndex(U, S);
	auto rbond = commonIndex(S, V);
	auto ldelta = delta(lbond, bb);
	auto rdelta = delta(rbond, prime(ba));
	S = copyTensor(S, IndexSet(bb, prime(ba)));
	U *= ldelta;
	V *= rdelta;
	lb = S/norm(S);
	auto inla = inverM(la);
	gb = U*inla;
	ga = V*inla;
}