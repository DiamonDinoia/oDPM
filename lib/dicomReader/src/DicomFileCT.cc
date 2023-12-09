//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include "DicomFileCT.hh"

std::vector<G4int> DicomFileCT::GetMaterislIDs() const{
    std::vector<G4int> mateIDs;
    G4int fCompress = theFileMgr->GetCompression();
    if( DicomFileMgr::verbose >= warningVerb ) G4cout << fLocation << " DumpMateIDsToTextFile "
                                                      << fFileName << " " << fMateIDs.size() << G4endl;
    for( int ir = 0; ir < fNoVoxelsY/fCompress; ir++ ) {
        for( int ic = 0; ic < fNoVoxelsX/fCompress; ic++ ) {
            mateIDs.emplace_back(fMateIDs[ic+ir*fNoVoxelsX/fCompress]);
        }
    }
    return mateIDs;
}
