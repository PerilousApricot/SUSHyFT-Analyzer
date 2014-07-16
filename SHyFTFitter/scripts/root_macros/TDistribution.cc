#include "TDistribution.h"


TDistribution::TDistribution()
{
  sample = "empty string";
  legentry = "empty string";
  filename = "empty string";
  var = "empty string";
  error = "";
  condition = "";

  line_color = 0;
  line_width = 0;
  fill_color = 0;
  marker_style = 0;
  marker_color = 0;
 }

TDistribution::~TDistribution()
{
}

