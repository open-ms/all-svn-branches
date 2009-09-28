#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <iostream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@brief a little theoretical ms/ms generator for some fixed given peptides.

	Generates ms/ms spectra for some fixed given peptides
	- no losses
	- no isotopes
	- fixed ion intensity 1
	- precursor (no peak, just precursor info in spectrum)
	- choosable is the mix of ion-types
	- choosable is the chargestate of the ms/ms spectra (preceeding (smaller) charge states are also added!)

	<B>The command line parameters of this tool are:</B>
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

	class TheoreticalMSMSGenerator
	: public TOPPBase
	{
	 public:
		TheoreticalMSMSGenerator()
			: TOPPBase("TheoreticalMSMSGenerator","little theoretical ms/ms generator for some fixed given peptides.")
		{

		}

 protected:

	typedef std::vector< FASTAFile::FASTAEntry > FASTAdata;

	virtual void registerOptionsAndFlags_()
	{
		registerOutputFile_("out","<file>","","output (simulated MS map) in mzML format",true);
		registerInputFile_("in","<file>","","Input protein sequences in FASTA format",true);
		registerFlag_("a","add a-ions");
		registerFlag_("b","add b-ions");
		registerFlag_("c","add c-ions");
		registerFlag_("x","add x-ions");
		registerFlag_("y","add y-ions");
		registerFlag_("z","add z-ions");
		registerIntOption_ ("charge", "Charge state for the ms/ms spectra", 2, "Charge state for the ms/ms spectra", false);
	}

	// Load proteins from FASTA file
	void loadFASTA_(const String& filename, FASTAdata& fastadata)
	{
		writeLog_(String("Loading sequence data from ") + filename +  String(" ..") );

		FASTAFile fastafile;

		// load FASTA file contents
		fastafile.load(filename, fastadata);
		writeLog_(String("done (") + fastadata.size() + String(" peptide(s) loaded)"));
	}

	ExitCodes main_(int, const char**)
	{
		//-------------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------
		String inputfile_name = getStringOption_("in");
		inputFileReadable_(inputfile_name);
		String outputfile_name = getStringOption_("out");
		outputFileWritable_(outputfile_name);
		bool aions(false),bions(false),cions(false),xions(false),yions(false),zions(false);
		if (getFlag_("a"))
		{
			aions = true;
		}
		if (getFlag_("b"))
		{
			bions = true;
		}
		if (getFlag_("c"))
		{
			cions = true;
		}
		if (getFlag_("x"))
		{
			xions = true;
		}
		if (getFlag_("y"))
		{
			yions = true;
		}
		if (getFlag_("z"))
		{
			zions = true;
		}

		Int charge(getIntOption_("charge"));


		//-------------------------------------------------------------
		// loading input
		//-------------------------------------------------------------
		FASTAdata fastadata;
		loadFASTA_(inputfile_name, fastadata);

		//-------------------------------------------------------------
		// generate a theoretical experiment map
		//-------------------------------------------------------------
		RichPeakMap my_rich_map;
		writeLog_("Starting generation");
		StopWatch w;

		w.start();

		//-------------------------------------------------------------
		// stuff
		//-------------------------------------------------------------
		TheoreticalSpectrumGenerator generator;
		Param p;
		bool losses = false; // "Neutral losses"
		p.setValue("add_losses", losses, "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
		bool isotopes = false; // "Isotope clusters"
		p.setValue("add_isotopes", isotopes, "If set to 1 isotope peaks of the product ion peaks are added");
		bool precursors = false;

		p.setValue("a_intensity", 1.0, "Intensity of the a-ions");
		p.setValue("b_intensity", 1.0, "Intensity of the b-ions");
		p.setValue("c_intensity", 1.0, "Intensity of the c-ions");
		p.setValue("x_intensity", 1.0, "Intensity of the x-ions");
		p.setValue("y_intensity", 1.0, "Intensity of the y-ions");
		p.setValue("z_intensity", 1.0, "Intensity of the z-ions");
		DoubleReal rel_loss_int = (DoubleReal)(10) / 100.0;
		p.setValue("relative_loss_intensity", rel_loss_int, "Intensity of loss ions, in relation to the intact ion intensity");
		generator.setParameters(p);

		/*-----#
			specs ms2 seqs hard-coded

			//mod
			sequences.push_back(String("(MOD:09998)ALDEK(MOD:00445)LFLI"));
			sequences.push_back(String("(MOD:09998)K(MOD:00445)LFSDISAI"));
			sequences.push_back(String("(MOD:09998)GLLPK(MOD:00445)SLYL"));
			sequences.push_back(String("(MOD:09998)IVIEAIHTV"));
			sequences.push_back(String("(MOD:09998)TLLDHIRTA"));
			sequences.push_back(String("(MOD:09999)ALDEK(MOD:00445)LFLI"));
			sequences.push_back(String("(MOD:09999)K(MOD:00445)LFSDISAI"));
			sequences.push_back(String("(MOD:09999)K(MOD:00445)LFTHDIML"));
			sequences.push_back(String("(MOD:09999)K(MOD:00445)LIIDREVV"));
			sequences.push_back(String("(MOD:09999)K(MOD:00445)VDDTFYYV"));
			sequences.push_back(String("(MOD:09999)GLLPK(MOD:00445)SLYL"));
			sequences.push_back(String("(MOD:09999)RVYEALYYV"));
			sequences.push_back(String("(MOD:09999)TLLDHIRTA"));

			//unmod
			sequences.push_back(String("ALDEKLFLI"));
			sequences.push_back(String("KLFSDISAI"));
			sequences.push_back(String("GLLPKSLYL"));
			sequences.push_back(String("IVIEAIHTV"));
			sequences.push_back(String("TLLDHIRTA"));
			sequences.push_back(String("ALDEKLFLI"));
			sequences.push_back(String("KLFSDISAI"));
			sequences.push_back(String("KLFTHDIML"));
			sequences.push_back(String("KLIIDREVV"));
			sequences.push_back(String("KVDDTFYYV"));
			sequences.push_back(String("GLLPKSLYL"));
			sequences.push_back(String("RVYEALYYV"));
			sequences.push_back(String("TLLDHIRTA"));

		 #-----*/

		/// @improvement catch no ion species chosen

		//read ms2 sequences of the experiment
		for(Size i = 0; i < fastadata.size(); ++i)
		{

			AASequence aa_sequence(fastadata[i].sequence);

			if (aa_sequence.isValid())
			{
				RichPeakSpectrum rich_spec;

				//~ for(Int ladung = 1; ladung <= charge; ++ladung)
				//~ {
				UInt ladung = charge;
					if (aions) // "A-ions"
					{
						generator.addPeaks(rich_spec, aa_sequence, Residue::AIon, ladung);
					}
					if (bions) // "B-ions"
					{
						generator.addPeaks(rich_spec, aa_sequence, Residue::BIon, ladung);
					}
					if (cions) // "C-ions"
					{
						generator.addPeaks(rich_spec, aa_sequence, Residue::CIon, ladung);
					}
					if (xions) // "X-ions"
					{
						generator.addPeaks(rich_spec, aa_sequence, Residue::XIon, ladung);
					}
					if (yions) // "Y-ions"
					{
						generator.addPeaks(rich_spec, aa_sequence, Residue::YIon, ladung);
					}
					if (zions) // "Z-ions"
					{
						generator.addPeaks(rich_spec, aa_sequence, Residue::ZIon, ladung);
					}
					if (precursors) // "Precursor"
					{
						generator.addPrecursorPeaks(rich_spec, aa_sequence, ladung);
					}
				//~ }

				Precursor prec;
				prec.setCharge(charge);
				prec.setMZ(aa_sequence.getMonoWeight(Residue::Full, charge)/double(charge)); //sum all aa internal weights plus h20 weight

					//~ cerr << aa_sequence.toString() << ": " << aa_sequence.getMonoWeight(Residue::Full, charge) << endl;
					//~ cerr << aa_sequence.toString() << ": " << aa_sequence.getMonoWeight(Residue::Internal, charge) << endl;

				prec.setIntensity(1);
				prec.setMetaValue("prec_to", aa_sequence.toString());
				rich_spec.getPrecursors().insert(rich_spec.getPrecursors().begin(), prec);

				rich_spec.setRT(-1);
				rich_spec.setMSLevel(2);
				rich_spec.setName(String(aa_sequence.toString() + "(theoretical)"));

				//push the new spec in the theoretical experiment
				my_rich_map.push_back(rich_spec);
			}
			else
			{
				//else peptide sequence is invalid! maybe inform the user?
				return INCOMPATIBLE_INPUT_DATA;
			}

		}
		w.stop();
		writeLog_(String("Generation took ") + String(w.getClockTime()) + String(" seconds"));

		//save new experiment as mzML
		writeLog_(String("Storing generated map in: ") + outputfile_name);
		MzMLFile().store(outputfile_name, my_rich_map);

		return EXECUTION_OK;
	}
};

int main( int argc, const char** argv )
{
	TheoreticalMSMSGenerator tool;
	return tool.main(argc,argv);
}
/// @endcond
