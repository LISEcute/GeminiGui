/**
 *!\brief weighted Monte Carlo
 *
 * Class CWeight is a base class that deals with a weighted 
 * Monte Carlo scheme. It is used to enhance the probabilty of IMF emission.
 * To compensate for this, each event is given a weight.
 * This weight should be used when histogramming events.
 *
 */

class CWeight
{
 protected:

  double fact;    //!< weighting factor
  int iWeight;  //!< ==0, no weighting
  double runningWeight; //!< running weight of event
  void findFactor(double Glight, double Gimf, double Gfission, double Ggamma);
 public:

  int  chooseChannel(double Glight, double Gimf, double Gfission, double Ggamma, 
                                    double xran);
  void setWeightIMF();
  double getWeightFactor();
};
