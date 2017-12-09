generic

package regge.io is

   function str (source : in Integer;
                 width  : in Integer := 0) return string;

   function str (source : in Real;
                 width  : in Integer := 10) return string;

   function read_command_arg (the_arg : Integer) return Integer;

   procedure read_lattice;

   num : Integer;  -- which one of the 8 datasets to read/check

   simp01 : Array (1 .. n_simp1_max) of Array1dIntg (1 .. 2);  -- vertices of each leg
                                                               -- used only in the io

end regge.io;
