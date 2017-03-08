// gemc headers
#include "ftm_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

static ftmConstants initializeFTMConstants(int runno)
{
    ftmConstants ftmc;
    
    // do not initialize at the beginning, only after the end of the first event,
    // with the proper run number coming from options or run table
    if(runno == -1) return ftmc;
    
    int isec,ilay,istr;
    
    
    ftmc.sigma_td_max = 0   ;  // // very small transverse diffusion (temporary)
    ftmc.w_i          = 25.0;  // ionization energy

    ftmc.rmin  =  70.43;         // inner radius of disks
    ftmc.rmax  = 143.66;         // outer radius of disks
    ftmc.pitch =  0.560;         // pitch of the strips
    ftmc.nstrips = (int) floor(2.*(ftmc.rmax+ftmc.rmin)/ftmc.pitch); // Number of strips

    return ftmc;
}


void ftm_HitProcess::initWithRunNumber(int runno)
{
    if(this->ftmcc.runNo != runno)
    {
        cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
        ftmcc = initializeFTMConstants(runno);
        ftmcc.runNo = runno;
    }
}

map<string, double> ftm_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();

	// FTM ID:
	// layer, type, sector, strip

	int layer  = 2*identity[0].id + identity[1].id - 2 ;
	int sector = identity[2].id;
	int strip  = identity[3].id;

	if(verbosity>4)
	{
		trueInfos tInfos(aHit);
		cout <<  log_msg << " layer: " << layer << "  sector: " << sector << "  Strip: " << strip
		<< " x=" << tInfos.x << " y=" << tInfos.y << " z=" << tInfos.z << endl;
	}
	dgtz["hitn"]   = hitn;
	dgtz["layer"]  = layer;
	dgtz["sector"] = sector;
	dgtz["strip"]  = strip;

	return dgtz;
}



vector<identifier> ftm_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	double x, y, z;
	G4ThreeVector  xyz = aStep->GetPostStepPoint()->GetPosition(); //< Global Coordinates of interaction
    G4ThreeVector Lxyz = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(xyz); //< Local Coordinates of interaction
    
    x = Lxyz.x()/mm;
	y = Lxyz.y()/mm;
	z = Lxyz.z()/mm;

	vector<identifier> yid = id;
	class ftm_strip ftms;
//	ftms.fill_infos(Detector);

	int layer  = 2*yid[0].id + yid[1].id - 2 ;

	//yid[3].id = ftms.FindStrip(layer-1, x, y, z);
	double depe = aStep->GetTotalEnergyDeposit();
	//cout << "resolMM " << layer << " " << x << " " << y << " " << z << " " << depe << " " << aStep->GetTrack()->GetTrackID() << endl;
	vector<double> multi_hit = ftms.FindStrip(layer-1, x, y, z, depe, Detector, ftmcc);

	int n_multi_hits = multi_hit.size()/2;

	// closest strip
	//yid[4].id = (int) multi_hit[0];
	yid[3].id = (int) multi_hit[0];

	yid[0].id_sharing = multi_hit[1];
	yid[1].id_sharing = multi_hit[1];
	yid[2].id_sharing = multi_hit[1];
	yid[3].id_sharing = multi_hit[1];
	// yid[4].id_sharing = multi_hit[1];

	// additional strip
	for(int h=1; h<n_multi_hits; h++)
	{
		for(int j=0; j<3; j++)
		{
			identifier this_id;
			this_id.name       = yid[j].name;
			this_id.rule       = yid[j].rule;
			this_id.id         = yid[j].id;
			this_id.time       = yid[j].time;
			this_id.TimeWindow = yid[j].TimeWindow;
			this_id.TrackId    = yid[j].TrackId;
			this_id.id_sharing = multi_hit[3];
			yid.push_back(this_id);
		}
		// last id is strip
		identifier this_id;
		this_id.name       = yid[3].name;
		this_id.rule       = yid[3].rule;
		this_id.id         = (int) multi_hit[2];
		this_id.time       = yid[3].time;
		this_id.TimeWindow = yid[3].TimeWindow;
		this_id.TrackId    = yid[3].TrackId;
		this_id.id_sharing = multi_hit[3];
		yid.push_back(this_id);
	}

	return yid;
}


map< string, vector <int> >  ftm_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;

	return MH;
}

// - charge: returns charge/time digitized information / step
map< int, vector <double> > ftm_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double ftm_HitProcess :: voltage(double charge, double time, double forTime)
{
	return 0.0;
}


// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> ftm_HitProcess :: electronicNoise()
{
	vector<MHit*> noiseHits;

	// first, identify the cells that would have electronic noise
	// then instantiate hit with energy E, time T, identifier IDF:
	//
	// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);

	// push to noiseHits collection:
	// noiseHits.push_back(thisNoiseHit)

	return noiseHits;
}

// this static function will be loaded first thing by the executable
ftmConstants ftm_HitProcess::ftmcc = initializeFTMConstants(-1);


