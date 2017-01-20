



using namespace std;

#ifndef Site_h
#define Site_h 1


class Site
{
	public:
    Site();
    Site(int, int);
    void calculateEnergy(void);
    void totalEnergy(void);
    void changeConf(void);
    ~Site();
    static int numberOfSites;
    static double T;
    static double TEnergy,totalMag;
    int spin;
    static int L;
	void setPtr(void);
    void setLattice(Site** lattice,int L1,double T1){theLattice=lattice;L=L1;T=T1;};
	private:


    static Site** theLattice;

    double energy,dE;
    int i,j;
	Site* left;
    Site* right;
    Site* up;
    Site* down;

};

#endif
