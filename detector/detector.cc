// G4 headers
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4IntersectionSolid.hh"
#include "G4NistManager.hh"
#include "G4Para.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4EllipticalTube.hh"
#include "G4Paraboloid.hh"
#include "G4Hype.hh"
#include "G4Sphere.hh"
#include "G4GenericTrap.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4PVReplica.hh"
#include "G4UnitsTable.hh"
#include "G4RotationMatrix.hh"
#include "G4TwoVector.hh"

// gemc headers
#include "detector.h"
#include "utils.h"

// mlibrary
#include "gstring.h"
using namespace gstring;

// C++ headers
#include <string>
#include <vector>
using namespace std;

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

detector::detector()
{
	SolidV    = NULL;
	LogicV    = NULL;
	PhysicalV = NULL;
}

int detector::create_solid(goptions gemcOpt, map<string, detector> *Map)
{
	int built = 0;
	if(SolidV) delete SolidV;
	
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Solid: >> ";
	double VERB    = gemcOpt.optMap["G4P_VERBOSITY"].arg ;
	string catch_v = gemcOpt.optMap["CATCH"].args;

	// cad, gdml obects will not be built
	if(description.find("gdmlParsed") != string::npos) return 0;
	if(description.find("cadImported") != string::npos) return 0;


	if(type.find("ReplicaOf") != string::npos)
	{
		if(VERB>4 || name.find(catch_v) != string::npos)
			cout << hd_msg << " " << name << " is a Replica. Solid Volume will not be built." << endl;
		return 0;
	}
	

	// ############################
	// CopyPlacement:
	// Point LogicV to the original
	// ############################
	if(type.find("CopyOf") != string::npos && type.find("CopyOf") == 0)
	{
		hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Copy: >> ";
		string original(type, 6, 190);
		
		// Look for original
		map<string, detector>::iterator it = (*Map).find(trimSpacesFromString(original));
		if(it == (*Map).end())
		{
			cout <<  hd_msg << " <" << original << "> not found. Exiting." << endl << endl;
			exit(0);
		}
		else
		{
			if(VERB>4 || name.find(catch_v) != string::npos)
			{
				cout << hd_msg << " " << name << " is a copy of <" << trimSpacesFromString(original) << ">. Pointing to its logical volume." << endl;
			}
			SetLogical(it->second.GetLogical());
		}
		built = 1;
	}
	
	
	// ################
	// Solid Operations
	// ################
	if(type.find("Operation:") != string::npos && type.find("Operation:") == 0)
	{
		hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Operation: >> ";
		
		bool translationFirst = false;
		
		// If Operation:~ it will perform the translation first
		size_t posTld = type.find("~");
		if( posTld != string::npos )
		{
			translationFirst = true;
			type.replace( posTld, 1, " " );
		}
		
		//////////////////////////////////
		bool absolutecoordinates = false;
		
		// If Operation:@ it will assume that position of second object is given in the common mother volume of both objects
		
		size_t pos_at = type.find("@");
		if( pos_at != string::npos )
		{
		    absolutecoordinates = true;
		    type.replace( pos_at, 1, " " );
		}
		///////////////////////////////////////////////
		
		
		size_t posp = 0;
		size_t posm = 0;
		size_t post = 0;
		size_t pos  = 0;
		// harcoded max size of 200 here
		string operation(type, 10, 190);
		posp = operation.find("+");
		posm = operation.find("-");
		post = operation.find("*");
		if     (posp != string::npos) pos = posp;
		else if(posm != string::npos) pos = posm;
		else if(post != string::npos) pos = post;
		if(!posp && !posm && !post)
		{
			cout << hd_msg << " Operation " << operation << " for " << name << " not recognized. Exiting." << endl;
			exit(0);
		}
		
		// Locating solids
		string solid1, solid2;
		string tsolid1, tsolid2;
		solid1.assign(operation, 0,     pos);
		solid2.assign(operation, pos+1, operation.size());
		tsolid1 = trimSpacesFromString(solid1);
		tsolid2 = trimSpacesFromString(solid2);
		
		// Locating second solid transformation
		map<string, detector>::iterator it1 = (*Map).find(tsolid1);
		map<string, detector>::iterator it2 = (*Map).find(tsolid2);
		if(it1 == (*Map).end())
		{
			cout <<  hd_msg << " " << tsolid1 << " Not found. Exiting." << endl << endl;
			exit(0);
		}
		if(it2 == (*Map).end())
		{
			cout <<  hd_msg << " " << tsolid2 << " Not found. Exiting." << endl << endl;
			exit(0);
		}
		
		// Define rotational and translational transformations then combine them
		G4RotationMatrix rotate    = it2->second.rot ;
		G4ThreeVector    translate = it2->second.pos;
		G4RotationMatrix invRot    =  rotate.invert() ;
		G4Transform3D    transf1( invRot, G4ThreeVector( 0, 0, 0 ) );
		G4Transform3D    transf2( G4RotationMatrix(), translate );
		G4Transform3D    transform = transf2 * transf1 ;
		
		if( absolutecoordinates && trimSpacesFromString(it1->second.mother) == trimSpacesFromString(it2->second.mother) )
		{
			//assume that second object position and rotation are given in absolute (mother) coordinates:
			
			G4RotationMatrix invrot1 = (it1->second.rot).inverse();
			G4RotationMatrix rotate2 = it2->second.rot;
			
			G4ThreeVector net_translation = it2->second.pos - it1->second.pos;
			
			// The net rotation should be the INVERSE of the rotation of object 2 relative to object 1.
			// rotate1/2 is the rotation of object 1/2 relative to their common mother.
			// If R is the rotation of 2 relative to 1, then R2 = R * R1 --> R = R2 R1^-1
			G4RotationMatrix invnet_rotation = (rotate2 * invrot1).invert();
			
			// In order to express the relative position of object 2 in the coordinate system of object one, we must rotate it
			// applying the same rotation as that used to position object 1, according to the GEANT4 framework:
			// I do not quite understand WHY this works, but through trial and error, I have
			// discovered that the combination of operations below is what works:
			net_translation *= it1->second.rot;
			transform = G4Transform3D( invnet_rotation, net_translation );
			//We don't want there to be any possibility to overwrite "transform" in this special case, so we force translationFirst to false here:
			translationFirst = false;
		}
		
		
		// If there was tilda in the operation string then the rotation and translation are switched
		// with respect to the default behaviour of G4UnionSolid with separate rotational matrix
		// and translatin vector
		if( translationFirst )
		{
			transform = transf1 * transf2 ;
		}
		
		if(posp != string::npos)
		{
			SolidV = new G4UnionSolid(       name, it1->second.GetSolid(), it2->second.GetSolid(), transform );
		}
		if(posm != string::npos)
		{
			SolidV = new G4SubtractionSolid( name, it1->second.GetSolid(), it2->second.GetSolid(), transform );
		}
		if(post != string::npos)
		{
			SolidV = new G4IntersectionSolid(name, it1->second.GetSolid(), it2->second.GetSolid(), transform );
		}
		
		if(VERB>4 || name.find(catch_v) != string::npos)
		{
			cout << hd_msg << " " << name << " is the  " << (pos==posp ? " sum " : " difference ") << " of " << tsolid1 << " and " << tsolid2 << endl;;
		}
		built = 1;
	}
	
	
	
	
	if(VERB>4 || name.find(catch_v) != string::npos)
	{
		cout << hd_msg << " " << name << " solid " << type << " built." << endl;
	}
	
	
	if(built==0)
	{
		cout << hd_msg << " " << name << " solid >" << type << "< not recognized. Exiting." << endl;
		exit(0);
	}
	return 1;
}


int detector::create_logical_volume(map<string, G4Material*> *MMats, goptions gemcOpt)
{
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Logical: >> ";
	double VERB    = gemcOpt.optMap["G4P_VERBOSITY"].arg ;
	string catch_v = gemcOpt.optMap["CATCH"].args;
	string defmat  = gemcOpt.optMap["DEFAULT_MATERIAL"].args;

	// cad, gdml obects will not be built
	if(description.find("gdmlParsed") != string::npos) return 0;
	if(description.find("cadImported") != string::npos) return 0;

	vector<aopt> changeMatOptions = gemcOpt.getArgs("SWITCH_MATERIALTO");
	for (unsigned int f = 0; f < changeMatOptions.size(); f++)
	{
		vector < string > oldNewMats = getStringVectorFromStringWithDelimiter(changeMatOptions[f].args, ",");
		if(oldNewMats.size() == 2)
		{
			// oldNewMats[0] = old
			// oldNewMats[1] = new
			if(material == trimSpacesFromString(oldNewMats[0]))
				material = trimSpacesFromString(oldNewMats[1]);
		}
	}

	// don't build the logical volumes for components or replicas
	if(material == "Component" || material == "OfReplica")
	{
		if(VERB>4 || name.find(catch_v) != string::npos)
			cout << hd_msg << " " << name << " is a Solid Component or a Replicant. Logical Volume will not be built." << endl;
		return 0;
	}
	
	// Check if Material Exists
	map<string, G4Material*>::iterator i = MMats->find(material);
	
	// if material is not defined, look in G4 table.
	if(i == MMats->end() && LogicV == 0)
	{
		G4NistManager* matman = G4NistManager::Instance();
		if(matman->FindOrBuildMaterial(material)) (*MMats)[material] = matman->FindOrBuildMaterial(material);
	}
	i = MMats->find(material);
	
	// if material is still not defined, use air
	if(i == MMats->end() && LogicV == 0)
	{
		if(defmat == "none")
		{
			cout << hd_msg << " Warning: material >" << material << "< is not defined for volume >" << name <<"<. Exiting" << endl;
			cout << hd_msg << " You can set the DEFAULT_MATERIAL flag to replace an undefined material. " << endl;
			exit(0);
		}
		else
		{
			material = defmat;
			if(MMats->find(material)== MMats->end())
			{
				cout << hd_msg << " Warning: " << defmat << " set with DEFAULT_MATERIAL is not found. Exiting" << endl;
				exit(0);
				
			}
		}
	}
	
	// Logical Volume Basic Constructor
	// If LogicV exists already, this is a copy
	if(LogicV == 0)
		LogicV = new G4LogicalVolume(SolidV, (*MMats)[material], name, 0, 0, 0, true);
	
	if(name == "root") LogicV->SetVisAttributes(G4VisAttributes::GetInvisible());
	else LogicV->SetVisAttributes(VAtts);
	
	if(VERB>4 || name.find(catch_v) != string::npos)
	{
		cout << hd_msg << " " << name << " Logical Volume built. " << endl;
	}


	return 1;
}


int detector::create_physical_volumes(goptions gemcOpt, G4LogicalVolume *mamma)
{
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Physical: >> ";
	double VERB    = gemcOpt.optMap["G4P_VERBOSITY"].arg ;
	bool   OVERL   = gemcOpt.optMap["CHECK_OVERLAPS"].arg > 0 ;
	string catch_v = gemcOpt.optMap["CATCH"].args;
	if(PhysicalV) delete PhysicalV;

	// cad, gdml obects will not be built
	if(description.find("gdmlParsed") != string::npos) return 0;
	if(description.find("cadImported") != string::npos) return 0;

	// don't build physical volumes for components or replicas.
	// Replicas are built in the dedicated routine
	if(material == "Component" || material == "OfReplica")
	{
		if(VERB>4 || name.find(catch_v) != string::npos)
			cout << hd_msg << " " << name << " is a Solid Component or a Replicant. Physical Volume will not be built, or replicas will be built instead." << endl;
		return 1;
	}
	
	
	if(name == "root")
		PhysicalV = new G4PVPlacement(0,          ///< rotation
									  G4ThreeVector(),   ///< translation
									  LogicV,            ///< logical volume
									  name.c_str(),      ///< name
									  0,                 ///< Mother Logical Volume
									  false,             ///< no boolean operation
									  0);                ///< copy number
	
	else
		PhysicalV = new G4PVPlacement(&rot,          ///< rotation
									  pos,                  ///< translation
									  LogicV,               ///< logical volume
									  name.c_str(),         ///< name
									  mamma,                ///< Mother Logical Volume
									  false,                ///< pMany (for future use)
									  ncopy,                ///< ncopy
									  OVERL);               ///< Checks Volume Overlapping at Placement time

	if(VERB>4 || name.find(catch_v) != string::npos)
	{
		if(mamma)
			cout << hd_msg << " " << name << " Physical Volume(s) built inside " << mamma->GetName() << "." << endl;
	}
	
	return 1;
}

int detector::create_replicas(goptions gemcOpt, G4LogicalVolume *mamma, detector replicant)
{
	string hd_msg  = gemcOpt.optMap["LOG_MSG"].args + " Physical: >> ";
	double VERB    = gemcOpt.optMap["G4P_VERBOSITY"].arg ;
	string catch_v = gemcOpt.optMap["CATCH"].args;
	
	// deleting the replicant detector physicsal volume as well
	if(replicant.PhysicalV) delete replicant.PhysicalV;
	
	// reset this physical volume
	if(PhysicalV) delete PhysicalV;
	
	
	if(type.find("ReplicaOf:") == 0)
	{
		if(dimensions.size() != 4)
		{
			cout << hd_msg << " Fatal Error: the number of parameters for " << name
				  << " is " << dimensions.size() <<  ":" << endl;
			for(unsigned int i=0; i<dimensions.size(); i++)
				cout << "      parameter " << i + 1 << ": " <<  dimensions[i] << endl;
			cout << "      This does not match a G4Replicas (4). Exiting" << endl << endl;
			exit(0);
		}
		
		EAxis pAxis;
		if(dimensions[0] == 1) pAxis = kXAxis;
		if(dimensions[0] == 2) pAxis = kYAxis;
		if(dimensions[0] == 3) pAxis = kZAxis;
		int nreps     = (int) get_number(dimensions[1]);
		double width  = get_number(dimensions[2]);
		double offset = get_number(dimensions[3]);
		
		// the logical volume is built
		// the mother is built as well
		PhysicalV = new G4PVReplica(name,   ///< name
											 replicant.LogicV, ///< Logical Volume to be replicated
											 mamma,            ///< Mother Logical Volume
											 pAxis,            ///< EAxis of copy - can be kXAxis,kYAxis,kZAxis,kRho,kRadial3D,kPhi
											 nreps,            ///< Number of repetitions
											 width,            ///< Width of repetitions
											 offset);          ///< Offset of repetitions
		
	}
	
	if(VERB>4 || name.find(catch_v) != string::npos)
	{
		if(mamma)
			cout << hd_msg << " " << name << " Physical Volume(s) built inside " << mamma->GetName() << "." << endl;
	}
	
	return 1;
}


ostream &operator<<(ostream &stream, detector Detector)
{
	cout  << endl;
	cout << "   Detector name:  "  << Detector.name        << "  -  " <<  Detector.description << endl;
	cout << "   Mother:  "         << Detector.mother                                 << endl;
	cout << "   Position (cm):  "  << Detector.pos/cm                                 << endl;
	cout << "   Rotation:       "  << Detector.rot                                    << endl;
	cout << "   Color:  "          << Detector.VAtts.GetColour()                      << endl;
	cout << "   Type:  "           << Detector.type                                   << endl;
	vector< vector<string> > dtypes = dimensionstype( Detector.type.c_str() );
	
	if(dtypes.size() != Detector.dimensions.size() && Detector.type.find("CopyOf") != 0)
	{
		for(unsigned int i=0; i<Detector.dimensions.size(); i++)
			cout << "   Size " << i + 1 << ":  "  << Detector.dimensions[i]  << endl;
	}
	if(dtypes.size() == Detector.dimensions.size() && Detector.type.find("CopyOf") != 0)
		for(unsigned int i=0; i<dtypes.size(); i++)
			cout << "   " << dtypes[i][0] << ":  "  << G4BestUnit(Detector.dimensions[i], dtypes[i][1])  << endl;
	
	cout << "   Material:  "       << Detector.material                               << endl;
	cout << "   Magnetic Field:  " << Detector.magfield                               << endl;
	cout << "   Copy Number:  "    << Detector.ncopy                                  << endl;
	cout << "   Activated: "       << ( Detector.exist==1 ?   "yes"   : "no" )        << endl;
	cout << "   Visible: "         << ( Detector.visible==1 ? "yes"   : "no" )        << endl;
	cout << "   Style: "           << ( Detector.style==1 ?   "solid" : "wireframe" ) << endl;
	cout << "   Sensitivity: "     << Detector.sensitivity                            << endl;
	if(Detector.sensitivity != "no")
		cout << "   hitType: "       << Detector.hitType                               << endl;
	
	if(Detector.identity.size())
		cout << Detector.identity ;

	if(Detector.SolidV) {
		cout << "   Volume: "       << bestValueUnits(Detector.SolidV->GetCubicVolume(), "Volume") << endl;
		cout << "   Surface Area: " << bestValueUnits(Detector.SolidV-> GetSurfaceArea(), "Surface") << endl;
	}
	if(Detector.LogicV) {
		cout << "   Mass: " << bestValueUnits(Detector.LogicV->GetMass(), "Mass") << endl;
	}
	cout << endl;
	
	return stream;
}

bool detector::operator == (const detector& D) const
{
	// Name uniquely identifies a volume
	if(D.name == this->name)
		return true;
	else return false;
}


























