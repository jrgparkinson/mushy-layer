/*
 * Channel.h
 *
 *  Created on: 28 Jul 2017
 *      Author: parkinsonjl
 */

#ifndef SRC_CHANNEL_H_
#define SRC_CHANNEL_H_

#include "IntVectSet.H"

/// Representation of a brine channel
/**
 * Consider a channel as just a collection of IntVects (an IntVectSet)
 * but with a few other helpful functions e.g. calculating the width
 */
class Channel : public IntVectSet
{
public:
  //  Channel ();

  /// Standard creator
  Channel();

  /// Create with intvectset
  Channel(const IntVect &iv) : IntVectSet(iv) { m_finished = false; }

  /// Does this channel border the intvect?
  bool borders(const IntVect &iv);

  /// Get the width at the vertical position specified
  Real width(Side::LoHiSide a_side, int a_offset, Real a_dx);

  /// Channel height
  Real height(Real a_dx);

  /// Chanel width, averaged over height
  Real averageWidth(Real a_dx);

  /// Horizontal position in grid
  /**
   * Implicit assumption here is all channels are symmetric
   */
  Real location();

  /// Remove bottom row of cells
  void removeBottomCells();

  /// Have we finished defining a channel
  bool isFinished();

  /// Set isFinished to true
  void setFinished();

  /// Destructor
  virtual ~Channel();

  /// Compute channel spacing
  static void channelSpacing(Vector<Real> &a_spacing,
                             Vector<Channel *> &a_channels, Real a_dx,
                             ProblemDomain a_probDomain)
  {
    int numChannels = a_channels.size();
    int numSpacings =
        a_probDomain.isPeriodic(0) ? numChannels : numChannels - 1;

    if (numSpacings <= 0)
    {
      a_spacing.resize(0);
      return;
    }
    else
    {
      a_spacing.resize(numChannels);
    }

    Vector<Real> channelLocations(numChannels);

    for (int i = 0; i < numChannels; i++)
    {
      channelLocations[i] = a_channels[i]->location();
    }

    channelLocations.sort();

    // Get all spacings that don't involve wrapping around a periodic domain
    for (int i = 0; i < (numChannels - 1); i++)
    {
      a_spacing[i] = (channelLocations[i + 1] - channelLocations[i]) * a_dx;
    }

    // If periodic, get the one extra channel spacing
    if (a_probDomain.isPeriodic(0))
    {
      int domMin_i = a_probDomain.domainBox().smallEnd(0);
      int domMax_i = a_probDomain.domainBox().bigEnd(0);

      a_spacing[numSpacings - 1] =
          ((domMax_i - channelLocations[numChannels - 1]) +
           (channelLocations[0] - domMin_i)) *
          a_dx;
    }
  }

private:
  /// Set if we've finished defining a channel
  bool m_finished;
};

#endif /* SRC_CHANNEL_H_ */
