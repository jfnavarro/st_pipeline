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

import org.apache.hadoop.io.RawComparator;
import org.apache.hadoop.io.WritableComparator;
import org.apache.hadoop.io.WritableUtils;


/**
 * Group based only on the sequence location.
 */
public class GroupByLocationComparator implements RawComparator<SequenceId>
{
	@Override
	public int compare(byte[] b1, int s1, int l1, byte[] b2, int s2, int l2) {

		int sizeVint1 = WritableUtils.decodeVIntSize(b1[s1]);
		int sizeVint2 = WritableUtils.decodeVIntSize(b2[s2]);

		int retval = WritableComparator.compareBytes(b1, s1+sizeVint1, l1-sizeVint1-Byte.SIZE/8, // Byte.SIZE is the byte size in bits
																								 b2, s2+sizeVint2, l2-sizeVint2-Byte.SIZE/8);
		return retval;
	}

	@Override
	public int compare(SequenceId s1, SequenceId s2) {
		return s1.getLocation().compareTo(s2.getLocation());
	}
}
