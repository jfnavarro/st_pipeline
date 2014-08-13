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

public enum AlignFlags {
	Paired(1),
	ProperlyPaired(2),
	Unmapped(4),
	MateUnmapped(8),
	OnReverse(16),
	MateOnReverse(32),
	Read1(64),
	Read2(128),
	SecondaryAlignment(256),
	FailedQC(0x200),
	Duplicate(0x400);

	private final int bit;

	private AlignFlags(int bit) {
		this.bit = bit;
	}

	public boolean is(int flag) {
		return (bit & flag) != 0;
	}

	public boolean isNot(int flag) {
		return (bit & flag) == 0;
	}

	/**
	 * Set this flag's bit in bitset.
	 */
	public int set(int bitset) {
		return bitset | bit;
	}

	/**
	 * Clear this flag's bit from bitset.
	 */
	public int clear(int bitset) {
		return bitset & ~bit;
	}
}
