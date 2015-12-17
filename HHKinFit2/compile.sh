 if [[ -z "$KINFIT2PATH" ]]
 then
     echo "Trying to initialize environment:"
     source ./setup.sh
 fi
 if [[ -z "$KINFIT2PATH" ]]
 then
     echo "Trying to initialize environment failed!"
 else
     echo "removing old files"
     rm -f libHHKinFit2.so
    
     echo "creating shared library"
     g++ -Wall -fPIC -shared src/*.cpp `root-config --cflags --glibs` -I $KINFIT2PATH -D HHKINFIT2 -o libHHKinFit2.so 
     echo "creating executable"
#     #g++ bin/main.C `root-config --cflags --glibs` -I $KINFIT2PATH -L $KINFIT2PATH/HHKinFit2/HHKinFit2 -lHHKinFit2  -o runHHKinFit
#     #g++ -std=c++11 bin/compareKinFits.C `root-config --cflags --glibs` -I $KINFIT2PATH -L $KINFIT2PATH/HHKinFit2/HHKinFit2 -lHHKinFit2 -L $KINFIT1PATH -lHHKinFit1 -o compareKinFits
#     g++ -std=c++11 bin/compare2DKinFits.C `root-config --cflags --glibs` -I $KINFIT2PATH -L $KINFIT2PATH/HHKinFit2/HHKinFit2 -lHHKinFit2 -L $KINFIT1PATH -lHHKinFit1 -o compare2DKinFits
#     #g++ -std=c++11 bin/gensimHH.C `root-config --cflags --glibs` -I $KINFIT2PATH -L $KINFIT2PATH/HHKinFit2/HHKinFit2 -lHHKinFit2 -o gensimHH
#     #g++ -std=c++11 bin/example.C `root-config --cflags --glibs` -I $KINFIT2PATH -L $KINFIT2PATH/HHKinFit2/HHKinFit2 -lHHKinFit2 -o example
    g++ bin/fitSingleEvent2DKinFit.C `root-config --cflags --glibs` -I $KINFIT2PATH -L $KINFIT2PATH/HHKinFit2/HHKinFit2 -lHHKinFit2 -D HHKINFIT2 -o fitSingleEvent2DKinFit
    g++ bin/fitSingleEvent1DKinFit.C `root-config --cflags --glibs` -I $KINFIT2PATH -L $KINFIT2PATH/HHKinFit2/HHKinFit2 -lHHKinFit2 -D HHKINFIT2 -o fitSingleEvent1DKinFit
 fi

