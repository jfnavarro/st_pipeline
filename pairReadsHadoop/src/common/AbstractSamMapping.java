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

import java.util.List;

/**
 * Abstract readable mapping that contains the data in a SAM record.
 *
 * This class provides interpretation of the SAM CIGAR string.  It does not
 * enforce any kind of storage format or medium.
 */
public abstract class AbstractSamMapping extends AbstractTaggedMapping
{
	////////////////////////////////////////////////
	// variables
	////////////////////////////////////////////////

	/**
	 * Cache the alignment returned by getAlignment().
	 */
	protected List<AlignOp> alignment;

	////////////////////////////////////////////////
	// methods
	////////////////////////////////////////////////

	/**
	 * Get the entire text related to tag so it can be scanned.
	 * For instance, given a SAM record with tags
	 * ...XT:A:U	NM:i:0	SM:i:37	AM:i:0	X0:i:1...
	 * getTagText("NM") should return "NM:i:0"
	 * while
	 * getTagText("XX") should return null.
	 * @return null if the tag isn't found.  The tag's text otherwise.
	 */
	abstract protected String getTagText(String name);

	/**
	 * Scan a SAM tag to produce an AbstractTaggedMapping.TagCacheItem.
	 */
	protected TagCacheItem makeTagItem(String name) throws NoSuchFieldException
	{
		String text = getTagText(name);
		if (text == null)
			throw new NoSuchFieldException("no tag with name " + name);
		String[] fields = text.split(":", 3);
		if (fields.length < 3)
			throw new FormatException("Invalid SAM tag syntax " + text);
		if (fields[1].length() != 1)
			throw new FormatException("Invalid SAM tag type syntax: " + text);

		return new TagCacheItem(TagDataType.fromSamType(fields[1].charAt(0)), fields[2]);
	}

	public List<AlignOp> getAlignment() throws IllegalStateException
	{
		if (isUnmapped())
			throw new IllegalStateException();

		// scan the CIGAR string and cache the results
		if (alignment == null)
			alignment = AlignOp.scanCigar(getCigarStr());

		return alignment;
	}

	abstract public String getCigarStr() throws IllegalStateException;
}
