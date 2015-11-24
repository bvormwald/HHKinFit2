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
     g++ -Wall -fPIC -shared src/*.cpp `root-config --cflags --glibs` -I ./interface -D HHKINFIT2 -o libHHKinFit2.so 
     echo "creating executable"
#     #g++ examples/main.C `root-config --cflags --glibs` -I $KINFIT2PATH/include -L $KINFIT2PATH -lHHKinFit2  -o runHHKinFit
#     #g++ -std=c++11 examples/compareKinFits.C `root-config --cflags --glibs` -I $KINFIT2PATH/include -I $KINFIT1PATH/interface -L $KINFIT2PATH -lHHKinFit2 -L $KINFIT1PATH -lHHKinFit1 -o compareKinFits
#     g++ -std=c++11 examples/compare2DKinFits.C `root-config --cflags --glibs` -I $KINFIT2PATH/include -I $KINFIT1PATH/interface -L $KINFIT2PATH -lHHKinFit2 -L $KINFIT1PATH -lHHKinFit1 -o compare2DKinFits
#     #g++ -std=c++11 examples/gensimHH.C `root-config --cflags --glibs` -I $KINFIT2PATH/include -I $KINFIT1PATH/interface -L $KINFIT2PATH -lHHKinFit2 -o gensimHH
#     #g++ -std=c++11 examples/example.C `root-config --cflags --glibs` -I $KINFIT2PATH/include -I $KINFIT1PATH/interface -L $KINFIT2PATH -lHHKinFit2 -o example
    g++ examples/fitSingleEvent2DKinFit.C `root-config --cflags --glibs` -I $KINFIT2PATH/interface -L $KINFIT2PATH -lHHKinFit2 -D HHKINFIT2 -o fitSingleEvent2DKinFit
    g++ examples/fitSingleEvent1DKinFit.C `root-config --cflags --glibs` -I $KINFIT2PATH/interface -L $KINFIT2PATH -lHHKinFit2 -D HHKINFIT2 -o fitSingleEvent1DKinFit
 fi

