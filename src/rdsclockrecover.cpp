void (std::vector<float> &rds, std::vector<int> &sym){

  //rds is the sampled sinusoid as input
  //syms is the high low and low high combination needed for manchester decoding

  float min = rds[0];
  int mini = 0;


  //initiliaze this based on sample rate to skip once local min/max found
  int skip = 0;

  for(int i=1; i <rds.size(); i++){

    //find local min / local max
    if(fabs(rds[i-1]-rds[i]) < min){

      min = fabs(rds[i-1]-rds[i];
      mini = i;

    }

    if(fabs(rds[mini])>fabs(rds[mini-1]) && fabs(rds[mini])>fabs(rds[mini+1])){

      if(rds[mini]>0){

        sym[0]  = 1;

      }

      else{

        sym[0]  = 0;

      }

      break;

    }
  }

  int count = mini + skip;

  for(int i = 1; count<rds.size(); i++){

    if(rds[count]>0){

      sym[i]  = 1;

    }

    else{

      sym[i]  = 0;

    }

    count += skip;

  }


}
