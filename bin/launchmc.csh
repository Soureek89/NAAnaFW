python new_singletop.py --t3batch -f files/pisanew/renamed/ -P SingleElectron -S 14 -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -d DATA -c electron
python new_singletop.py --t3batch -f files/pisanew/renamed/ -P SingleElectron -S 14 -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -d DATA -c electronantiiso

python new_singletop.py --t3batch -f files/pisanew/renamed/ -P SingleMuon -S 14 -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -d DATA -c muon
python new_singletop.py --t3batch -f files/pisanew/renamed/ -P SingleMuon -S 14 -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -d DATA -c muonantiiso


#MC

python new_singletop.py --t3batch -f files/pisanew/renamed/ -P ST,VJ -S 5 -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -c muon,electron
python new_singletop.py --t3batch -f files/pisanew/renamed/ -P TT -S 3 -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -c muon,electron
python new_singletop.py --t3batch -f files/pisanew/renamed/ -P VV -S 5 -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -c muon,electron

python new_singletop.py --t3batch -f files/pisanew/renamed/ -P ST,VJ -S 5 -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -c muonantiiso,electronantiiso
python new_singletop.py --t3batch -f files/pisanew/renamed/ -P TT -S 3 -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -c muonantiiso,electronantiiso
python new_singletop.py --t3batch -f files/pisanew/renamed/ -P VV -S 5 -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -c muonantiiso,electronantiiso

#python new_singletop.py --t3batch -f files/pisanew/renamed/ -P VV,QCDMu -S 7 -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -c muon
#python new_singletop.py --t3batch -f files/pisanew/renamed/ -P VV,QCDMu -S 7 -t $WORKSPACE/public/xWajid/trees -o $WORKSPACE/public/xWajid/ -c muonantiiso







