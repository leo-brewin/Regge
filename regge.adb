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

with Ada.Numerics;                              use Ada.Numerics;
with Ada.Numerics.Generic_Elementary_Functions;

package body regge is

   package Maths is new Ada.Numerics.Generic_Elementary_Functions (Real); use Maths;

   function "*" (Left : Integer; Right : Real) return Real is
   begin
      return Real(Left) * Right;
   end "*";

   function "*" (Left : Real; Right : Integer) return Real is
   begin
      return Left * Real(Right);
   end "*";

   function sign (x : Real) return Integer is
   begin
      if x = 0.0 then
         return 0;
      else
         if x < 0.0
            then return -1;
            else return +1;
         end if;
      end if;
   end sign;

   function inv_angle
     (n1_dot_m2 : Real;
      m1_dot_m2 : Real;
      sig_m1    : Integer) return Real
   is
      phi : Real;
      rho : Real;
   begin

      case bone_signature is

         when timelike =>

            phi := arccos (m1_dot_m2);

         when spacelike =>

            if abs (m1_dot_m2) < abs (n1_dot_m2)
               then rho := sign (n1_dot_m2) * m1_dot_m2;
               else rho := sign (m1_dot_m2) * n1_dot_m2;
            end if;

            phi := sig_m1 * arcsinh (rho);

      end case;

      return phi;

   end inv_angle;

   procedure gauss_system
     (solution   : out Array2dReal;
      the_matrix : in  Array2dReal;
      the_rhs    : in  Array2dReal;
      num_row    : in  Integer;
      num_rhs    : in  Integer;
      cond       : out Real;
      failed     : out Boolean)
   is

      singular_pivot : constant := 1.0e-20;

      -- Solves a set of linear equations by Gaussian elimination. ------------

      det                                    : Real;
      pivoting                               : Boolean;
      size                                   : Integer := num_row;
      col, row, pivot_row                    : Integer;
      sum, pivot, factor, maximum, temporary : Real;
      min_pivot, max_pivot                   : Real;
      rhs                                    : Array2dReal := the_rhs;
      matrix                                 : Array2dReal := the_matrix;

   begin

      col       := 1;
      row       := 1;
      det       := 1.0e0;
      min_pivot := 1.0e20;
      max_pivot := 0.0e0;

      failed   := False;
      pivoting := True;

      while pivoting loop

         -- search for pivot row

         maximum   := abs (matrix (col, col));
         pivot_row := col;
         for row in col + 1 .. size loop
            if maximum < abs (matrix (row, col)) then
               maximum   := abs (matrix (row, col));
               pivot_row := row;
            end if;
         end loop;

         -- swap diagonal row with pivot row

         if pivot_row /= col then
            for i in col .. size loop
               temporary             := matrix (pivot_row, i);
               matrix (pivot_row, i) := matrix (col, i);
               matrix (col, i)       := temporary;
            end loop;
            for i in 1 .. num_rhs loop
               temporary          := rhs (pivot_row, i);
               rhs (pivot_row, i) := rhs (col, i);
               rhs (col, i)       := temporary;
            end loop;
            det := -det;
         end if;

         -- set diagonal element to one

         pivot := matrix (col, col);

         if abs (pivot) < min_pivot then
            min_pivot := abs (pivot);
         end if;
         if abs (pivot) > max_pivot then
            max_pivot := abs (pivot);
         end if;

         if abs (pivot / max_pivot) < singular_pivot then
            pivoting := False;
            failed   := True;
            det      := 0.0e0;

         else
            matrix (col, col) := 1.0e0;
            for i in col + 1 .. size loop
               matrix (col, i) := matrix (col, i) / pivot;
            end loop;
            for i in 1 .. num_rhs loop
               rhs (col, i) := rhs (col, i) / pivot;
            end loop;

            -- eliminate elements below diagonal

            for row in col + 1 .. size loop
               factor := matrix (row, col);
               for i in col .. size loop
                  matrix (row, i) := matrix (row, i) - factor * matrix (col, i);
               end loop;
               for i in 1 .. num_rhs loop
                  rhs (row, i) := rhs (row, i) - factor * rhs (col, i);
               end loop;
            end loop;
            col      := col + 1;
            pivoting := (not failed) and (col <= size);
            det      := det * pivot;
         end if;

      end loop;

      -- form the back-substitution

      if not failed then
         for i in 1 .. num_rhs loop
            solution (size, i) := rhs (size, i);
            for row in reverse 1 .. size - 1 loop
               sum := 0.0e0;
               for col in row + 1 .. size loop
                  sum := sum + matrix (row, col) * solution (col, i);
               end loop;
               solution (row, i) := rhs (row, i) - sum;
            end loop;
         end loop;
      else
         for i in 1 .. size loop
            for j in 1 .. num_rhs loop
               solution (i, j) := 0.0e0;
            end loop;
         end loop;
      end if;

      if not failed
         then cond := min_pivot / max_pivot;
         else cond := 0.0e0;
      end if;

   end gauss_system;

   function inverse
     (matrix : Array2dReal) return Array2dReal
   is
      cond : Real;
      failed : Boolean;
      sol : Array2dReal (matrix'first(1)..matrix'last(1),matrix'first(1)..matrix'last(1));
      identity : Array2dReal (matrix'first(1)..matrix'last(1),matrix'first(1)..matrix'last(1));
   begin
      if matrix'length(1) /= matrix'length (2)
         then raise Constraint_Error with " matrix must be square";
      else
         for i in identity'range(1) loop
            identity (i,i) := 1.0;
            for j in i+1 .. identity'last(2) loop
               identity (i,j) := 0.0;
               identity (j,i) := 0.0;
            end loop;
         end loop;
         gauss_system (sol,matrix,identity,matrix'length(1),matrix'length(2),cond,failed);
         if not failed
            then return sol;
            else raise Constraint_Error with " matrix may be singular";
         end if;
      end if;
   end inverse;

   procedure set_normal
     (n     : out Array1dReal;
      m     : out Array1dReal;
      sig_n : out Integer;
      sig_m : out Integer;
      c     :     Integer;
      d     :     Integer;
      inv_metric : Array2dReal)
   is
      tmp  : Real;
   begin

      -- compute the outward pointing unit normal, n

      if inv_metric (c, c) > 0.0 then

         sig_n := +1;
         tmp   := +Sqrt (abs (inv_metric (c, c)));
         for a in 1 .. 4 loop
            n (a) := -inv_metric (a, c) / tmp;  -- minus sign for outward normal
         end loop;

      else

         sig_n := -1;
         tmp   := -Sqrt (abs (inv_metric (c, c)));
         for a in 1 .. 4 loop
            n (a) := -inv_metric (a, c) / tmp;  -- minus sign for outward normal
         end loop;

      end if;

      -- compute the unit tangent vector, m

      tmp := inv_metric (c, d) / inv_metric (c, c);

      for a in 1 .. 4 loop
         m (a) := inv_metric (a, d) - tmp * inv_metric (a, c);
      end loop;

      if m (d) > 0.0 then

         sig_m := +1;
         tmp   := +Sqrt (abs (m (d)));
         for a in 1 .. 4 loop
            m (a) := m (a) / tmp;
         end loop;

      else

         sig_m := -1;
         tmp   := -Sqrt (abs (m (d)));
         for a in 1 .. 4 loop
            m (a) := m (a) / tmp;
         end loop;

      end if;

   end set_normal;

   function get_signature
     (bone : Integer) return bone_type
   is
      leg_ij, leg_ik, leg_jk : Integer;
      lsq_ij, lsq_ik, lsq_jk : Real;
      g11, g22, g12 : Real;
   begin

      leg_ij := simp12 (bone)(1);
      leg_ik := simp12 (bone)(2);
      leg_jk := simp12 (bone)(3);

      lsq_ij := lsq (leg_ij);
      lsq_ik := lsq (leg_ik);
      lsq_jk := lsq (leg_jk);

      g11 := lsq_ij;
      g22 := lsq_ik;
      g12 := (g11 + g22 - lsq_jk) / 2.0;

      if g11*g22-g12*g12 > 0.0
         then return spacelike;
         else return timelike;
      end if;

   end get_signature;

   procedure set_metric
     (metric : out Array2dReal;
      bone   :     Integer;
      index  :     Integer)
   is

      offset : Integer;

      leg_ab : Integer;
      leg_ij, leg_ik, leg_jk : Integer;
      leg_ai, leg_aj, leg_ak : Integer;
      leg_bi, leg_bj, leg_bk : Integer;

      lsq_ab : Real;
      lsq_ij, lsq_ik, lsq_jk : Real;
      lsq_ia, lsq_ja, lsq_ka : Real;
      lsq_ib, lsq_jb, lsq_kb : Real;

   begin

      offset := 4*index - 4;

      leg_ij := simp12 (bone)(1);
      leg_ik := simp12 (bone)(2);
      leg_jk := simp12 (bone)(3);

      leg_ai := simp12 (bone)(offset+4);
      leg_aj := simp12 (bone)(offset+5);
      leg_ak := simp12 (bone)(offset+6);

      leg_ab := simp12 (bone)(offset+7);

      leg_bi := simp12 (bone)(offset+8);
      leg_bj := simp12 (bone)(offset+9);
      leg_bk := simp12 (bone)(offset+10);

      lsq_ij := lsq (leg_ij);
      lsq_ik := lsq (leg_ik);
      lsq_jk := lsq (leg_jk);

      lsq_ia := lsq (leg_ai);
      lsq_ja := lsq (leg_aj);
      lsq_ka := lsq (leg_ak);

      lsq_ib := lsq (leg_bi);
      lsq_jb := lsq (leg_bj);
      lsq_kb := lsq (leg_bk);

      lsq_ab := lsq (leg_ab);

      metric (1, 1) := lsq_ij;
      metric (2, 2) := lsq_ik;
      metric (3, 3) := lsq_ia;
      metric (4, 4) := lsq_ib;

      metric (1, 2) := (metric (1, 1) + metric (2, 2) - lsq_jk) / 2.0;
      metric (1, 3) := (metric (1, 1) + metric (3, 3) - lsq_ja) / 2.0;
      metric (1, 4) := (metric (1, 1) + metric (4, 4) - lsq_jb) / 2.0;
      metric (2, 3) := (metric (2, 2) + metric (3, 3) - lsq_ka) / 2.0;
      metric (2, 4) := (metric (2, 2) + metric (4, 4) - lsq_kb) / 2.0;
      metric (3, 4) := (metric (3, 3) + metric (4, 4) - lsq_ab) / 2.0;

      metric (2, 1) := metric (1, 2);
      metric (3, 1) := metric (1, 3);
      metric (4, 1) := metric (1, 4);
      metric (3, 2) := metric (2, 3);
      metric (4, 2) := metric (2, 4);
      metric (4, 3) := metric (3, 4);

   end set_metric;

   procedure get_defect
     (defect   : out Real;
      bone     :     Integer;
      n_vertex :     Integer)
   is

      phi : Real;

      procedure take_one_step
        (phi   : in out Real;
         index :        Integer)
      is
         n1_dot_m2  : Real;
         m1_dot_m2  : Real;
         tmp_a      : Real;
         tmp_b      : Real;
         tmp_c      : Real;
         sig_m1     : Integer;
         metric     : Array2dReal (1..4,1..4);
         inv_metric : Array2dReal (1..4,1..4);
      begin

         set_metric (metric, bone, index);                        -- standard metric for (i,j,k,a,b);
         inv_metric := inverse (metric);

         -- update phi, the defect

         tmp_a := inv_metric (3, 3) * inv_metric (4, 4);
         tmp_b := inv_metric (3, 4) * inv_metric (3, 4);
         tmp_c := 1.0 - (tmp_b / tmp_a);

         sig_m1 := sign (tmp_c * inv_metric (3, 3));

         n1_dot_m2 := -sign (inv_metric (4, 4)) * Sqrt (abs (tmp_c));
         m1_dot_m2 := -sign (tmp_a - tmp_b) * inv_metric (3, 4) / Sqrt (abs (tmp_a));

         phi := phi + inv_angle (n1_dot_m2, m1_dot_m2, sig_m1);

      end take_one_step;

   begin

      phi := 0.0;

      -- compute the defect

      -- take a tour around the bone, visiting each 4-simplex in turn

      for a in 1 .. n_vertex loop
         take_one_step (phi, a);
      end loop;

      case bone_signature is
         when timelike  => defect := 2.0 * Pi - phi;
         when spacelike => defect := -phi;
      end case;

   end get_defect;

   procedure get_defect_deriv
     (deriv    : out Real;
      bone     :     Integer;
      n_vertex :     Integer;
      leg_pq   :     Integer)
   is

      rho : Real;
      deriv_metric : Array2dReal (1 .. 4, 1 .. 4);

      procedure set_metric_deriv
        (bone  : Integer;
         index : Integer)
      is
         the_delta  : Real;
         save_lsq   : Real;
         metric_p   : Array2dReal (1 .. 4, 1 .. 4);
         metric_m   : Array2dReal (1 .. 4, 1 .. 4);
      begin

         -- compute derivative of the metric
         -- since the metric is linear in the lsq, a simple finite difference
         -- method will give the exact answer

         save_lsq  := lsq (leg_pq);
         the_delta := save_lsq / 10.0;

         lsq (leg_pq) := save_lsq + the_delta;   set_metric (metric_p, bone, index);
         lsq (leg_pq) := save_lsq - the_delta;   set_metric (metric_m, bone, index);

         for a in 1 .. 4 loop
            for b in 1 .. 4 loop
               deriv_metric (a, b) := (metric_p (a, b) - metric_m (a, b)) / (2.0 * the_delta);
            end loop;
         end loop;

         -- reset lsq to its original value

         lsq (leg_pq) := save_lsq;

      end set_metric_deriv;

      procedure take_one_step
        (rho   : in out Real;
         index :        Integer)
      is
         sum            : Real;
         sig_n1, sig_m1 : Integer;
         sig_n2, sig_m2 : Integer;
         n1, m1         : Array1dReal (1 .. 4);
         n2, m2         : Array1dReal (1 .. 4);
         metric         : Array2dReal (1..4,1..4);
         inv_metric     : Array2dReal (1..4,1..4);
      begin

         set_metric (metric, bone, index);                        -- standard metric for (i,j,k,a,b);
         inv_metric := inverse (metric);

         set_normal (n1, m1, sig_n1, sig_m1, 4, 3, inv_metric);   -- outward normal & tangent vector to (i,j,k,a)
         set_normal (n2, m2, sig_n2, sig_m2, 3, 4, inv_metric);   -- outward normal & tangent vector to (i,j,k,b)

         set_metric_deriv (bone, index);

         -- update rho, the derivative of the defect

         sum := 0.0;
         for a in 1 .. 4 loop
            for b in 1 .. 4 loop
               sum := sum + (n1 (a) * m1 (b) + n2 (a) * m2 (b)) * deriv_metric (a, b);
            end loop;
         end loop;

         rho := rho + sum;

      end take_one_step;

   begin

      rho := 0.0;

      for a in 1 .. n_vertex  loop
         take_one_step (rho, a);
      end loop;

      case bone_signature is
         when timelike  => deriv := -rho / 2.0;
         when spacelike => deriv :=  rho / 2.0;
      end case;

   end get_defect_deriv;

   procedure get_defect_and_deriv
     (defect   : out Real;
      deriv    : out Real;
      bone     :     Integer;
      n_vertex :     Integer;
      leg_pq   :     Integer)
   is

      phi : Real;
      rho : Real;
      deriv_metric : Array2dReal (1 .. 4, 1 .. 4);

      procedure set_metric_deriv
        (bone  : Integer;
         index : Integer)
      is
         the_delta  : Real;
         save_lsq   : Real;
         metric_p   : Array2dReal (1 .. 4, 1 .. 4);
         metric_m   : Array2dReal (1 .. 4, 1 .. 4);
      begin

         -- compute derivative of the metric
         -- since the metric is linear in the lsq, a simple finite difference
         -- method will give the exact answer

         save_lsq  := lsq (leg_pq);
         the_delta := save_lsq / 10.0;

         lsq (leg_pq) := save_lsq + the_delta;   set_metric (metric_p, bone, index);
         lsq (leg_pq) := save_lsq - the_delta;   set_metric (metric_m, bone, index);

         for a in 1 .. 4 loop
            for b in 1 .. 4 loop
               deriv_metric (a, b) := (metric_p (a, b) - metric_m (a, b)) / (2.0 * the_delta);
            end loop;
         end loop;

         -- reset lsq to its original value

         lsq (leg_pq) := save_lsq;

      end set_metric_deriv;

      procedure take_one_step
        (phi   : in out Real;
         rho   : in out Real;
         index :        Integer)
      is
         sum            : Real;
         n1_dot_m2      : Real;
         m1_dot_m2      : Real;
         tmp_a          : Real;
         tmp_b          : Real;
         tmp_c          : Real;
         sig_n1, sig_m1 : Integer;
         sig_n2, sig_m2 : Integer;
         n1, m1         : Array1dReal (1 .. 4);
         n2, m2         : Array1dReal (1 .. 4);
         metric         : Array2dReal (1..4,1..4);
         inv_metric     : Array2dReal (1..4,1..4);
      begin

         set_metric (metric, bone, index);                        -- standard metric for (i,j,k,a,b);
         inv_metric := inverse (metric);

         set_normal (n1, m1, sig_n1, sig_m1, 4, 3, inv_metric);   -- outward normal & tangent vector to (i,j,k,a)
         set_normal (n2, m2, sig_n2, sig_m2, 3, 4, inv_metric);   -- outward normal & tangent vector to (i,j,k,b)

         set_metric_deriv (bone, index);

         -- update phi, the defect

         tmp_a := inv_metric (3, 3) * inv_metric (4, 4);
         tmp_b := inv_metric (3, 4) * inv_metric (3, 4);
         tmp_c := 1.0 - (tmp_b / tmp_a);

         sig_m1 := sign (tmp_c * inv_metric (3, 3));

         n1_dot_m2 := -sign (inv_metric (4, 4)) * Sqrt (abs (tmp_c));
         m1_dot_m2 := -sign (tmp_a - tmp_b) * inv_metric (3, 4) / Sqrt (abs (tmp_a));

         phi := phi + inv_angle (n1_dot_m2, m1_dot_m2, sig_m1);

         -- update rho, the derivative of the defect

         sum := 0.0;
         for a in 1 .. 4 loop
            for b in 1 .. 4 loop
               sum := sum + (n1 (a) * m1 (b) + n2 (a) * m2 (b)) * deriv_metric (a, b);
            end loop;
         end loop;

         rho := rho + sum;

      end take_one_step;

   begin

      phi := 0.0;
      rho := 0.0;

      for a in 1 .. n_vertex  loop
         take_one_step (phi, rho, a);
      end loop;

      case bone_signature is
         when timelike  => defect := 2.0 * Pi - phi;
         when spacelike => defect := -phi;
      end case;

      case bone_signature is
         when timelike  => deriv := -rho / 2.0;
         when spacelike => deriv :=  rho / 2.0;
      end case;

   end get_defect_and_deriv;

end regge;
