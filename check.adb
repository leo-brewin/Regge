with Ada.Text_IO;   use Ada.Text_IO;
with regge;
with regge.io;

procedure check is

	type Real is digits 18;

   package regge_lib    is new regge (Real);   use regge_lib;
   package regge_lib_io is new regge_lib.io;   use regge_lib_io;

   package Real_IO    is new Ada.Text_IO.Float_IO (Real);      use Real_IO;
   package Integer_IO is new Ada.Text_IO.Integer_IO (Integer); use Integer_IO;

   txt    : File_Type;

   defect : Real;
   deriv  : Real;
   bone   : Integer;

begin

   num := read_command_arg (1);

   read_lattice;

   bone := 1;

   bone_signature := get_signature (bone);

   Create (txt,Out_file,"results/data-0"&str(num,1)&".txt");

   get_defect (defect, bone, n_loop02(bone));

   put_line (txt,"# defect");
   put_line (txt, str(defect,25));
   new_line (txt);
   put_line (txt,"# i  j   d(defect)/d(lsq(i,j))");

   for a in 1..n_simp1 loop

      get_defect_deriv (deriv, bone, n_loop02(bone), a);

      put_line (txt, str (simp01(a)(1),3) &' '&
                     str (simp01(a)(2),2) &' '&
                     str (deriv,25));

   end loop;

   Close (txt);

end check;
