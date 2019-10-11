	List<Pair<word> > mergePatchPairs = List<Pair<word> > (2);
	mergePatchPairs[0].first() = "region0/Right";
	mergePatchPairs[0].second() = "region1/Left";

	Info<< "Creating merge patch pairs" << nl << endl;

	// Create and add point and face zones and mesh modifiers
	List<pointZone*> pz(mergePatchPairs.size());
	List<faceZone*> fz(3*mergePatchPairs.size());
	List<cellZone*> cz(0);

//	forAll(mergePatchPairs, pairI)
//	{
		const word mergeName
                (
                    mergePatchPairs[pairI].first()
                  + mergePatchPairs[pairI].second()
                  + name(pairI)
                );

                pz[pairI] = new pointZone
                (
                    mergeName + "CutPointZone",
                    labelList(0),
                    0,
                    multiMesh[myRank]->mesh->pointZones()
                );

                // Master patch
                const word masterPatchName(mergePatchPairs[pairI].first());
                const polyPatch& masterPatch =
                    multiMesh[myRank]->mesh->boundaryMesh()[masterPatchName];

                labelList isf(masterPatch.size());

                forAll(isf, i)
                {
                    isf[i] = masterPatch.start() + i;
                }

                fz[3*pairI] = new faceZone
                (
                    mergeName + "MasterZone",
                    isf,
                    boolList(masterPatch.size(), false),
                    0,
                    multiMesh[myRank]->mesh->faceZones()
                );

                // Slave patch
                const word slavePatchName(mergePatchPairs[pairI].second());
                const polyPatch& slavePatch =
                    multiMesh[myRank]->mesh->boundaryMesh()[slavePatchName];

                labelList osf(slavePatch.size());

                forAll(osf, i)
                {
                    osf[i] = slavePatch.start() + i;
                }

                fz[3*pairI + 1] = new faceZone
                (
                    mergeName + "SlaveZone",
                    osf,
                    boolList(slavePatch.size(), false),
                    1,
                    multiMesh->mesh->faceZones()
                );

                // Add empty zone for cut faces
                fz[3*pairI + 2] = new faceZone
                (
                    mergeName + "CutFaceZone",
                    labelList(0),
                    boolList(0, false),
                    2,
                    multiMesh[myRank]->mesh->faceZones()
                );
//	}  // end of all merge pairs

        Info<< "Adding point and face zones" << endl;
        multiMesh[myRank]->mesh->addZones(pz, fz, cz);


        Info<< "Creating attachPolyTopoChanger" << endl;
        attachPolyTopoChanger polyMeshAttacher(*mesh);
        polyMeshAttacher.setSize(mergePatchPairs.size());

//	forAll(mergePatchPairs, pairI)
//	{
                const word mergeName
                (
                    mergePatchPairs[pairI].first()
                  + mergePatchPairs[pairI].second()
                  + name(pairI)
                );

                // Add the sliding interface mesh modifier
                polyMeshAttacher.set
                (
                    pairI,
                    new slidingInterface
                    (
                        "couple" + name(pairI),
                        pairI,
                        polyMeshAttacher,
                        mergeName + "MasterZone",
                        mergeName + "SlaveZone",
                        mergeName + "CutPointZone",
                        mergeName + "CutFaceZone",
                        mergePatchPairs[pairI].first(),
                        mergePatchPairs[pairI].second(),
                        slidingInterface::INTEGRAL, // always integral
                        intersection::VISIBLE
                    )
                );
//	}

	polyMeshAttacher.attach(true);
