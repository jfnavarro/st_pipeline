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

import java.util.ArrayList;
import java.util.List;
import java.util.regex.*;

public class AlignOp {

	public static enum Type {
		Match,
		Insert,
		Delete,
		SoftClip,
		HardClip,
		Skip,
		Pad;

		public char getSymbol()
		{
			switch (this) {
			case Match:
				return 'M';
			case Insert:
				return 'I';
			case Delete:
				return 'D';
			case SoftClip:
				return 'S';
			case HardClip:
				return 'H';
			case Skip:
				return 'N';
			case Pad:
				return 'P';
			default:
				throw new RuntimeException("BUG!  Missing cigar operator in AlignOp::Type::getSymbol.");
			}
		}

		public static Type fromSymbol(String sym) {
			if (sym.length() != 1)
				throw new IllegalArgumentException("Unrecognized alignment operation symbol " + sym);
			else
				return fromSymbol(sym.charAt(0));
		}

		public static Type fromSymbol(char sym) {
			switch (sym) {
			case 'M':
				return Type.Match;
			case 'I':
				return Type.Insert;
			case 'D':
				return Type.Delete;
			case 'S':
				return Type.SoftClip;
			case 'H':
				return Type.HardClip;
			case 'N':
				return Type.Skip;
			case 'P':
				return Type.Pad;
			default:
				throw new IllegalArgumentException("Unrecognized alignment operation symbol " + sym);
			}
		}
	}

	protected static final Pattern CigarElementPattern = Pattern.compile("(\\d+)([MIDNSHP])");

	private Type op;
	private int len;

	public AlignOp(Type op, int len) {
		this.op = op;
		this.len = len;
	}

	public AlignOp.Type getType() { return op; }
	public int getLen() { return len; }

	public boolean equals(Object other)
	{
		if (other instanceof AlignOp)
		{
			AlignOp otherAlignment = (AlignOp) other;
			return otherAlignment.op == this.op && otherAlignment.len == this.len;
		}
		else
			return false;
	}

	public String toString() { return "(" + op + "," + len + ")"; }

	/**
	 * Compute the SAM-style CIGAR string for the given alignment.
	 */
	public static String cigarStr(List<AlignOp> alignment)
	{
		if (alignment == null || alignment.isEmpty())
			return "*";

		StringBuilder builder = new StringBuilder(50);

		for (AlignOp op: alignment)
			builder.append(op.getLen()).append(op.getType().getSymbol());

		return builder.toString();
	}

	public static List<AlignOp> scanCigar(String cigar)
	{
		if ("*".equals(cigar))
			return new ArrayList<AlignOp>(0);
		else
		{
			ArrayList<AlignOp> result = new ArrayList<AlignOp>(5);
			Matcher m = CigarElementPattern.matcher(cigar);

			int lastPositionMatched = 0;
			while (m.find())
			{
				result.add( new AlignOp(AlignOp.Type.fromSymbol(m.group(2)), Integer.parseInt(m.group(1))) );
				lastPositionMatched = m.end();
			}

			if (lastPositionMatched < cigar.length())
				throw new FormatException("Invalid CIGAR pattern " + cigar);

			if (result.isEmpty())
				throw new FormatException("Unable to parse any alignments from CIGAR pattern " + cigar);

			return result;
		}
	}
}
