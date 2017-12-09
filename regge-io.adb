with Ada.Text_IO;                               use Ada.Text_IO;
with Ada.Strings;                               use Ada.Strings;
with Ada.Strings.Fixed;                         use Ada.Strings.Fixed;
with Ada.Command_Line;                          use Ada.Command_Line;

package body regge.io is

   package Real_IO    is new Ada.Text_IO.Float_IO (Real);                      use Real_IO;
   package Integer_IO is new Ada.Text_IO.Integer_IO (Integer);                 use Integer_IO;

   function str (source : in Integer;
                 width  : in Integer := 0) return string is
      result : string(1..width);
      wide_result : string(1..Integer'width);
   begin
      if width = 0 then
         Put(wide_result,source);
         return trim(wide_result,left);  -- flush left, returns a string just large enough to contain the integer
      else
         Put(result,source);
      return result;
      end if;
   end str;

   function str (source : in Real;
                 width  : in Integer := 10) return string is
      result : string(1..width);
   begin
      -- 4932 = largest exponent for 18 dec. digits
      -- so may need up to 4 digits in the exponent + 1 for the sign = 5
      if source = 0.0 then
         Put (result,source,width-7,3);
      elsif abs (source) < 1.0 then
         if abs (source) >= 1.0e-99 then
            Put (result,source,width-7,3);
         elsif abs (source) >= 1.0e-999 then
            Put (result,source,width-8,4);
         else
            Put (result,source,width-9,5);
         end if;
      else
         if abs (source) < 1.0e100 then
            Put (result,source,width-7,3);
         elsif abs (source) < 1.0e1000 then
            Put (result,source,width-8,4);
         else
            Put (result,source,width-9,5);
         end if;
      end if;
      return result;
   end str;

   function read_command_arg (the_arg : Integer) return Integer is
      last : Integer;
      the_arg_value : Integer;
   begin
      get (Ada.Command_Line.Argument (the_arg),the_arg_value,last);
      return the_arg;
   end read_command_arg;

   procedure read_lattice
   is
      txt   : File_Type;
      tmp   : Integer;
      bone  : Integer;
   begin

      Open (txt,In_file,"data/lattice/simp12.txt");

      Get (txt,n_simp2);   skip_line (txt);
      for a in 1..n_simp2 loop
         Get (txt,bone);             skip_line (txt);
         Get (txt,n_loop02(bone));   skip_line (txt);
         Get (txt,n_simp12(bone));   skip_line (txt);
         for b in 1..n_simp12(bone) loop
            Get (txt,tmp);
            Get (txt,simp12(bone)(tmp));
         end loop;
      end loop;

      Close (txt);

      Open (txt,In_file,"data/lsq/data-0"&str(num,1)&".txt");

      Get (txt,n_simp1);
      skip_line (txt);  -- skip trailing text on first line
      skip_line (txt);  -- skip second line
      for i in 1 .. n_simp1 loop
         Get (txt, simp01(i)(1));
         Get (txt, simp01(i)(2));
         Get (txt, lsq(i));
      end loop;

      Close (txt);

   end read_lattice;

end regge.io;
