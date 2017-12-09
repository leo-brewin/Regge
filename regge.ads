-- Copyright (c) 2017, Leo Brewin <Leo.Brewin@monash.edu>
--
-- Permission to use, copy, modify, and/or distribute this software for any
-- purpose with or without fee is hereby granted, provided that the above
-- copyright notice and this permission notice appear in all copies.
--
-- THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
-- WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
-- MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
-- ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
-- WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
-- ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
-- OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

generic

   type Float is digits <>;

   n_simp1_max  : Integer := 64;  -- max. number of legs in the lattice
   n_simp2_max  : Integer := 1;   -- max. number of bones in the lattice
   n_loop02_max : Integer := 10;  -- max. number of vertices in the loop around a bone

package regge is

   subtype Real is Float;

   type bone_type is (timelike, spacelike);

   type Array1dReal is array (Integer range <>) of Real;
   type Array2dReal is array (Integer range <>, Integer range <>) of Real;

   type Array1dIntg is array (Integer range <>) of Integer;

	function "*" (Left : Integer; Right : Real) return Real;
	function "*" (Left : Real; Right : Integer) return Real;

   function get_signature
     (bone : Integer) return bone_type;

   procedure set_metric
     (metric : out Array2dReal;
      bone   :     Integer;
      index  :     Integer);

   procedure get_defect
     (defect   : out Real;          -- the defect
      bone     :     Integer;       -- the bone
      n_vertex :     Integer);      -- the number of vertices that enclose the bone

   procedure get_defect_deriv
     (deriv    : out Real;          -- the derivative of the defect
      bone     :     Integer;       -- the bone
      n_vertex :     Integer;       -- the number of vertices that enclose the bone
      leg_pq   :     Integer);      -- use this leg to compute the derivative

   procedure get_defect_and_deriv
     (defect   : out Real;          -- the defect
      deriv    : out Real;          -- the derivative of the defect
      bone     :     Integer;       -- the bone
      n_vertex :     Integer;       -- the number of vertices that enclose the bone
      leg_pq   :     Integer);      -- use this leg to compute the derivative

   bone_signature : bone_type;

   n_simp12_max : Integer := 3 + 4*n_loop02_max;   -- max. number of legs in the loop around a bone

   n_lsq_max    : Integer := n_simp1_max;          -- max. number of legs in the lattice

   n_simp1      : Integer := 0;                    -- number of legs in the lattice
   n_simp2      : Integer := 0;                    -- number of bones in the lattice

   n_loop02     : Array1dIntg (1 .. n_simp2_max);  -- number of vertices in the loop around a bone
   n_simp12     : Array1dIntg (1 .. n_simp2_max);  -- number of legs in the loop around a bone

   lsq          : Array1dReal (1..n_lsq_max);      -- the lsq's of the lattice

   simp12       : Array (1 .. n_simp2_max) of Array1dIntg (1 .. n_simp12_max);  -- all legs for each bone

end regge;
