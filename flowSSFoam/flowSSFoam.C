#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    Info<< "\nStarting simulation\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        laminarTransport.correct();
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End of simulation\n" << endl;

    return 0;
}