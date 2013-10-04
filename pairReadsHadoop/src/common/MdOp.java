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

public class MdOp
{
	protected static final Pattern MatchPattern = Pattern.compile("\\d++");
	protected static final Pattern MismatchPattern = Pattern.compile("[AGCTN]++");
	protected static final Pattern DeletePattern = Pattern.compile("\\^[AGCTN]++");

	public static enum Type {
		Match,
		Mismatch,
		Delete;
	}

	private Type op;
	private int len;
	private String seq; // Empty string for match.
	                    // For mismatch, bases different from reference and for
											// deletion the bases deleted from reference.
	                    // note: does not include '^' character

	public MdOp(MdOp.Type op, int len) {
		this.op = op;
		this.len = len;
		this.seq = "";
	}

	public MdOp(MdOp.Type op, int len, String seq) {
		if (op == Type.Match && !seq.isEmpty())
			throw new IllegalArgumentException("non-empty sequence " + seq + " provided for Match operator");
		this.op = op;
		this.len = len;
		this.seq = seq;
	}

	public MdOp.Type getType() { return op; }
	public int getLen() { return len; }
	public String getSeq() { return seq; }

	public boolean equals(Object other)
	{
		if (other instanceof MdOp)
		{
			MdOp otherMd = (MdOp) other;
			return otherMd.op == this.op && otherMd.len == this.len &&
				((this.seq == null && otherMd.seq == null) ||
				 (this.seq != null && otherMd.seq != null && otherMd.seq.equals(this.seq)));
		}
		else
			return false;
	}

	public String toString() { return "(" + op + "," + len + "," + seq + ")"; }

	/**
	 * Scan an MD tag into a list of MdOp elements.
	 */
	public static List<MdOp> scanMdTag(String tag) throws FormatException
	{
		ArrayList<MdOp> result = new ArrayList<MdOp>(5);
		int length;
		int end = tag.length();

		Matcher m = MatchPattern.matcher(tag);
		if (!m.lookingAt())
			throw new FormatException("Invalid MD tag '" + tag + "'. Tag doesn't start with a number.");

		length = Integer.parseInt(m.group());

		if (length > 0)
			result.add(new MdOp(Type.Match, length));
		// else don't add a 0-length op

		// advance the scanner
		m.region(m.end(), end);

		while (tag.length() > m.regionStart())
		{
			m.usePattern(MismatchPattern);
			if (m.lookingAt()) // found a mismatch
			{
				result.add(new MdOp(Type.Mismatch, m.group().length(), m.group()));
				// advance the scanner
				m.region(m.end(), end);
			}
			else
			{
				m.usePattern(DeletePattern);
				if (m.lookingAt()) // found a deletion
				{
					result.add(new MdOp(Type.Delete, m.group().length() - 1, m.group().substring(1))); // -1 for the ^ character
					// advance the scanner
					m.region(m.end(), end);
				}
				else
					throw new FormatException("Invalid MD tag '" + tag + "' (pos " + m.regionStart() + "). Match number not followed by a mismatch or delete.");
			}

			m.usePattern(MatchPattern);
			if (m.lookingAt())
			{
				length = Integer.parseInt(m.group());
				if (length > 0)
					result.add(new MdOp(Type.Match, length));

				// advance the scanner
				m.region(m.end(), end);
			}
			else
				throw new FormatException("Invalid MD tag '" + tag + "' (pos " + m.regionStart() + "). Mismatch or delete not followed by a match number.");
		}
		return result;
	}
}
