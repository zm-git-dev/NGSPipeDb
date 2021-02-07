# Basic knowledge of Linux
1.The absolute path is identified by a forward slash (/), which is different from the Windows operating system<br />
2.How to use linux command: enter the command at the shell prompt: "$" . Command like: [Path]/command [-option parameter] [file|directory]<br />
  For example: ~/miniconda3/envs/RNA/bin/fastp -i /home/sysbio/SRR8467686.fastq.gz<br />
3.The cd command is used to switch the directory where the session is located. Command like: cd directory<br />
  For example: cd /usr/bin<br />
4.The pwd command is used to view the directory where the current session is located, command like: pwd<br />
  For example: pwd<br />
5.Enter a single dot (.) to indicate the current directory, and enter a double dot (..) to indicate the parent directory of the current directory<br />
6.Use the man command to view the operation manual of all required commands, command like: man command name<br />
  For example: man cp<br />
7.Filter:<br />
  The filter refers to the name of the specified file and is used after many commands. Add a filter after the ls command, then only the information of the file will be displayed<br />
For example: ls -l my_file<br />
This command specifies the relevant information of the output file my_file<br />
   <Font color = red> When writing a filter, you can use a question mark (?) to represent a character, an asterisk (*) to represent zero or more characters, and use the tab key (Tab) to quickly complete the file name or directory First name</Font><br />
8.The which command is used to find the path of the command. Command like: which grep<br />
  For example: which fastp<br />

# About permission
## View the file permission
When you use **ls -l** command, the permissions will be marked at the beginning of the file, like: drwxrwxrwx. The first letter represents the type of file: d means directory, - means file, and l means soft link. The remaining 9 letters are grouped into 3, representing the permissions for the owner, the permissions for the owner's group, and the permissions for other groups. r means read-only; w means allowing user to modify the content of the file; x means allowing user to execute the file as a program; - means user having no such permission<br />
## Modify permissions: chmod command
<Font color = red>You must have the authority to operate files</Font><br />
Command like: chmod (user permissions) (group permissions) (other permissions) file<br />
Permission is expressed in octal code: r=4, w=2, x=1<br />
For example: chmod 755 test.txt<br />
The result is -rwxr-xr-x test.txt<br />

# Command about directory
1.Check what files are in the directory-list function<br />
(1) Basic list function: ls command<br />
  command like: ls various parameters directory name<br />
  --F distinguish between files and directories<br />
  --R recursive option to list files in subdirectories contained in the current directory<br />
  --l produces a long list of output results, including information about each file<br />
  --d lists only the information of the directory itself, not its contents<br />
  --i The inode number of the file, which is the unique identifier of the file<br />
  --a Display both hidden files and ordinary files<br /><Font color = red>The parameters can be combined and written, for example: ls -Fd</Font><br />
(2) The way to output the tree list: tree tool<br />
   The tool needs to be downloaded by yourself. command like: tree filter<br />
2.Create a directory: mkdir command<br />
  command like: mkdir directory name<br />
   To create multiple directories and subdirectories, you can use the -p parameter, for example: mkdir -p New_file/work/file1<br />
   <Font color = red>The -p parameter in the mkdir command can create missing parent directories as needed</Font><br />
3.Delete directory: rmdir command<br />
  command like: $ rmdir parameter directory name<br />
  --No parameters after rmdir delete empty directories<br />
  --rf delete all contents in the directory<br />
  --i Ask a question before deleting<br />
  <Font color = red>Because Linux does not have a recycle bin, you must add the -i parameter when deleting to confirm whether the deleted content is correct</Font><br />
# Command about file
1.Create a file: touch command
  Command like: touch parameter file name<br />
  -Create a new file without adding parameters after touch, if the file already exists, change the modification time<br />
  -a Change the access time of an existing file<br />
2.Delete files: rm command<br />
  Command like: rm -i file name<br />
  <Font color = red>Because Linux does not have a recycle bin, you must add the -i parameter when deleting to confirm whether the deleted content is correct</Font><br />
3.Copy files: cp command
  Command like: cp [parameter] [file name]<br />
  --i Force to ask if you need to overwrite existing files<br />
  --R Copy the entire directory recursively<br />
<Font color = red>To avoid errors, it is recommended to add the -i parameter</Font><br />
4.Rename and move files: mv command
(1) Rename<br />
  command like: mv [original file name] [new file name]<br />
(2) Move<br />
  command like: mv file name in the original path, target new path<br />
(3) Rename and move operations can be performed at the same time<br />
  For example: mv /home/picture/book /home/file/cook Move the book file in the picture folder to the file folder and rename it to cook<br />
5.View files
<1> View file type: file command<br />
 command like: file file name<br />
 You can check the file type, character encoding method, whether the file can be run, etc.<br />
<2> View the entire file<br />
    (1) cat command: display all data in the text file<br />
      command like: cat parameter file name<br />
   --No parameters after cat display content<br />
   --n display content after adding line number<br />
   --b displays only after adding line numbers to text content<br />
   --T does not display tabs, replace tabs with ^T to display<br />
    (2) More command: display the content of the text file, but it will stop after each page is displayed<br />
    (3) less command: display the content of text files and support more advanced interaction<br />
<3> View some files<br />
    (1) head command: view the beginning of the file<br />
      command like: head -n file name<br />
      n is the number of rows displayed<br />
    (2) tail command: view the end of the file<br />
      command like: tail -n file name<br />
      n is the number of rows displayed<br />
<4> Upload files, download files<br />
    (1) rz file name: upload file<br />
    (2) sz file name: download file<br />
# Command about process<br />
1.Probe the process<br />
(1) ps command<br />
command like: ps parameter<br />
Symbol description:<br />
    -PID process ID<br />
    -Terminal device when TTY process starts<br />
    -TIME Cumulative CPU time required by the process<br />
    -The name of the program started by CMD<br />
(2) top command<br />
command like: top<br />
Symbol description:<br />
    -load average The three values are the average load of the last 1min, 5min, and 15min. The larger the value, the larger the load, and the more than 2 indicates the system is busy.<br />
    -PR process priority<br />
    -Moderate value of NI process<br />
    -The total amount of virtual memory occupied by the VIRT process<br />
    -The total amount of physical memory occupied by the RES process<br />
    -%MEM The ratio of memory used by the process to available memory<br />
    -S process state (D: interruptible sleep state, R: running, S: sleeping, T: tracking or stopping state, Z: rigid state)<br />
(3) The difference between ps and top command<br />
    The ps command displays information at a specific time point, and the top command displays real-time information.<br />
2.End the process<br />
(1) kill command<br />
command like: kill PID or kill -s process signal<br />
(2) killall command<br />
command like: killall process name<br />
<Font color = red>Use wildcards carefully in the killall command</Font><br />
# Execute tasks to the background
Method 1: [Path]/command [-option parameter] [file|directory] &<br />
Method 2: nohup [path]/command [-option parameter] [file|directory] &<br />
Method 3: screen command<br />
Command like :<br />
(1)Create a screen: screen -dmS screen_test<br />
(2)View the screen: screen -list<br />
(3)Connect to the screen: screen -r screen_test<br />
Common task management commands:<br />
(1)jobs: View tasks, return task number n and process number<br />
(2)bg %n: Transfer task number n to the background<br />
(3)fg %n: Transfer task number n to the foreground<br />
(4)ctrl+z: Suspend the current task<br />
(5)ctrl+c: Stop the current task<br />
(6)kill -n task: End the task<br />
# Command about disk space<br />
1.Mount the device: mount command<br />
command like: mount parameter file device type device to be mounted target location<br />
File device types are:<br />
-vfat: Windows long file system<br />
-ntfs: an advanced file system widely used in Windows<br />
-iso9660: Standard CD-ROM file system<br />
2.Uninstall the device: umount command<br />
command like: umount location/target device<br />
3.Explore the disk situation<br />
(1) df command<br />
command like: df [parameter] [target disk]<br />
(2) du command<br />
command like: du [parameter] [target disk]<br />
   --c displays the total size of all listed files<br />
   --h Display in a user-readable way<br />
   --s displays the total of each output parameter<br />
(3) The difference between df and du commands<br />
The df command displays the disk usage, the du command displays the disk usage of each file<br />
4.Disk partition: fdisk command
command like: fdisk -l [disk name]<br />
5.Disk formatting: mkfs command
command like: mksf -t file [system format] [parameter] [disk name]<br />
File system formats are: mkfs.cramfs, mkfs.ext2, mkfs.ext3, mkfs.msdos, mkfs.vfat<br />
6.Disk verification: fsck command
command like: fsck -t file [system format] [parameter] [disk name]<br />
File system formats are: fsck.cramfs, fsck.ext2, fsck.ext3, fsck.msdos, fsck.vfat<br />

# Case
1.Download software:<br />
wget –c ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz<br />
2.Unzip file: tar zxvf ncbi-blast-2.6.0+-x64-linux.tar.gz<br />
3.Install the software:<br />
(1)cd soft<br />
(2)./configure<br />
(3)make<br />
(4)make install<br />
4.Update PATH: <br />
(1) Add the ./ncbi-blast-2.6.0+/bin directory to the environment variables: vim ~/.bashrc<br />
(2)Add the following statement, save and exit: export PATH=./ncbi-blast-2.6.0+/bin:$PATH<br />
(3)Update environment variables: source ~/.bashrc<br />



```python

```
