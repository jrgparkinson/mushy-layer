/*
 * Channel.cpp
 *
 *  Created on: 28 Jul 2017
 *      Author: parkinsonjl
 */

#include "Channel.h"

Channel::Channel ()
{
  m_finished = false;
}

Channel::~Channel ()
{
}


bool Channel::borders(const IntVect& iv)
{
  Channel grownChannel(*this);
  grownChannel.grow(1);

  return grownChannel.contains(iv);
}

Real Channel::width(Side::LoHiSide a_side, int a_offset, Real a_dx)
{
  //  Box b = thisChannel.minBox();
  Box b;
  if (a_side == Side::Lo)
  {
    b = adjCellLo(minBox(), 1, 1);
    b.shift(1, 1);
  }
  else
  {
    b = adjCellHi(minBox(), 1, 1);
    b.shift(1, -1);
  }
  b.shift(1, a_offset);
  //          bottomBox.shift(1, 1);

  IntVectSet ivs(*this);
  ivs &= b;

  Real width = ivs.numPts()*a_dx;

  return width;
}

void Channel::removeBottomCells()
{
  Box bottomBox = adjCellLo(minBox(), 1, 1);
  bottomBox.shift(1, 1);
  minus_box(bottomBox);
}

Real Channel::height(Real a_dx)
{
  int cellHeight = minBox().bigEnd(1) - minBox().smallEnd(1);
  return cellHeight * a_dx;
}

Real Channel::averageWidth(Real a_dx)
{
  int maxOffset = minBox().bigEnd(1) + 1 - minBox().smallEnd(1);
  Real avWidth = 0.0;
  Side::LoHiSide startingSide = Side::Hi;
  for (int i = 0; i < maxOffset; i++)
  {
    avWidth += width(startingSide, -i, a_dx);
  }

  avWidth /= maxOffset;

  return avWidth;

}

Real Channel::location()
{
  // Calculate average x coordinate for all intvects
  Real average_i = 0;
  for (IVSIterator ivs(*this); ivs.ok(); ++ivs)
  {
    IntVect iv = ivs();
    average_i += iv[0];

  }

  average_i /= numPts();

  return average_i;
}

bool Channel::isFinished()
{
  return m_finished;
}

void Channel::setFinished()
{
  m_finished = true;
}


//void channelSpacing(Vector<Real>& a_spacing, Vector<Channel*>& a_channels, Real a_dx)

