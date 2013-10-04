// Copyright (C) 2011-2012 CRS4.
//
// This file is part of Seal.
//
// Seal is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// Seal is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License along
// with Seal.  If not, see <http://www.gnu.org/licenses/>.

package common;

import common.SequenceId;

import org.apache.hadoop.mapreduce.Partitioner;

/**
 * Partition based only on the sequence location.
 */
public class SequenceIdLocationPartitioner<V> extends Partitioner<SequenceId,V>
{
	@Override
	public int getPartition(SequenceId key, V value, int numPartitions)
	{
		// clear the sign bit with & Integer.MAX_VALUE instead of calling Math.abs,
		// which will return a negative number for Math.abs(Integer.MIN_VALUE).
		return (key.getLocation().hashCode() & Integer.MAX_VALUE) % numPartitions;
	}
}
