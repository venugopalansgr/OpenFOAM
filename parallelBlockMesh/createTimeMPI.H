//
// createTimeMPI.H
// ~~~~~~~~~~~~

    if (myRank==0) Foam::Info<< "Create time\n" << Foam::endl;

    Foam::Time runTime(Foam::Time::controlDictName, args);
