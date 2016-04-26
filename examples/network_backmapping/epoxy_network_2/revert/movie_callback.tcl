proc moviecallback { args } {
 if {$::MovieMaker::userframe > 200 && $::MovieMaker::userframe <= 2200} {
    mol modstyle 1 0 VDW [expr 2.0-(0.001*($::MovieMaker::userframe-200))] 62.000
 }
 animate next
}

proc setupmovie {} {
 enablemoviecallback
}

proc finishmovie {} {
  disablemoviecallback
}


## Easy-to-use proc to enable the user-defined movie frame callback
proc enablemoviecallback { }  {
  trace add variable ::MovieMaker::userframe write moviecallback
}

## Easy-to-use proc to disable the user-defined movie frame callback
proc disablemoviecallback { } {
  trace remove variable ::MovieMaker::userframe write moviecallback
}
