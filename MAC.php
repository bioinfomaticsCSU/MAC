#!/usr/bin/env php
<?php

function exchange(&$left, &$right){
	$temp = $left;
	$left = $right;
	$right = $temp;
}

function connect($left, $right, $coor1, $coor2){
	$temp1 = array_search($left[0], $right[$coor1]);
	$temp2 = array_search($left[1], $right[$coor2]);
	$arr1_temp = cycle($right[$coor1], $temp1);
	$arr2_temp = cycle($right[$coor2], $temp2);
	$arrConnect = array_merge($arr1_temp, $arr2_temp);
	unset($right[$coor1]);
	unset($right[$coor2]);
	$right[$coor1] = $arrConnect;
	return $right;
}

function cycle($right, $k){
	return array_merge(array_slice($right, $k), array_slice($right, 0, $k));
}

function multiply($left, $right){
	$left_mapping = array();
	$right_mapping = array();
	$temp = array();
	foreach($left as $i => $markers){
		$temp = $markers;
		array_shift($temp);
		array_push($temp, $markers[0]);
		$left_mapping += array_combine($markers, $temp);
	}
	$numMarkers = 0;
	$temp = array();
	foreach($right as $i => $markers){
		$temp = $markers;
		array_shift($temp);
		array_push($temp,$markers[0]);
		foreach($temp as $j => $n){
			$right_mapping[$markers[$j]] = $left_mapping[$n];
			$numMarkers++;
		}
	}
	$rod_temp = array();
	$temp = array();
	$i = 0;
	while($i < $numMarkers){
		foreach($right_mapping as $k => $marker){
			$m = isset($nextMarker) ? $nextMarker : $k;
			array_push($temp, $m);
			$i++;
			$nextMarker = $right_mapping[$m];
			unset($right_mapping[$m]);
			if(!isset($right_mapping[$nextMarker]) && $temp != array()){
				array_push($rod_temp, $temp);
				$temp = array();
				unset($nextMarker);
				break;
			}
		}
	}
	return $rod_temp;
}

function merge($array){
	$merge = array();
	foreach($array as $i){
		$merge = array_merge($merge, $i);
	}
	return $merge;
}

function parse_contig($multiFasta){
	$seqs = array();
	$parse = fopen($multiFasta, 'r');
	global $contigSeq;
	while(($seq = fgets($parse)) !== false){
		array_push($seqs, rtrim($seq));
	}
	fclose($parse);
	$pos = '';
	$neg = '';
	$forward_arr = array();
	$reverse_arr = array();
	$result = 0;
	foreach($seqs as $k => $seq){
		if($seq == '') continue;
		if($seq[0] != '>'){
			array_push($forward_arr, $seq."\n");
			array_push($reverse_arr, strrev($seq)."\n");
		}
		if($seq[0] == '>' || !isset($seqs[$k+1])){
			if($result == 1){
				foreach($forward_arr as $each){
					$pos .= $each;
				}
				$n = count($reverse_arr);
				for($i = $n-1; $i>=0; $i--){
					$neg .= $reverse_arr[$i];
				}
				$contigSeq[$fileName] = $pos;		
				$new_reverse = $neg;
				$seqLen = strlen($neg);
				for($i = 0; $i < $seqLen; $i++){
					if($neg[$i] == 'A'){
						$new_reverse[$i] = 'T';
					}
					else if($neg[$i] == 'T'){
						$new_reverse[$i] = 'A';
					}
					else if($neg[$i] == 'C'){
						$new_reverse[$i] = 'G';
					}
					else if($neg[$i] == 'G'){
						$new_reverse[$i] = 'C';
					}
					else if($neg[$i] == 'a'){
						$new_reverse[$i] = 't';
					}
					else if($neg[$i] == 't'){
						$new_reverse[$i] = 'a';
					}
					else if($neg[$i] == 'c'){
						$new_reverse[$i] = 'g';
					}
					else if($neg[$i] == 'g'){
						$new_reverse[$i] = 'c';
					}
				}
				$contigSeq["$fileName.reverse"] = $new_reverse;
				$forward_arr = array();
				$reverse_arr = array();
				$pos = '';
				$neg = '';
			}
			if(strstr($seq, ' ')){
				$fileName = substr($seq, 1, strpos($seq, ' ')-1);
			}
			else{
				$fileName = substr($seq, 1);
			}
			$result = 1;
		}
	}
}

function vertex_parse(&$vertex, $a, $b){
	$vex_new = $vertex[$b];
	$vals = array_keys($vertex, $vex_new);
	foreach($vals as $val){
		$vertex[$val] = $vertex[$a];
	}
}

function fastaFile($draftFileName, $outputFile){
	$arr = file($outputFile);
	global $contigSeq;
	global $plus_o;
	$att = '';
	$edge = 0;
	$tail = '.fasta';
	file_put_contents($outputFile.$tail, '');
	foreach($arr as $k => $each){
		if($each == "\n" || strstr($each, '>')){
			file_put_contents($outputFile.$tail, $each, FILE_APPEND);
			$edge = 1;
			continue;
		}
		if($edge == 0){
			$att = str_repeat('N', 100)."\n";
			file_put_contents($outputFile.$tail, $att, FILE_APPEND);
		}
		$partEach = trim(substr($each, 0, strpos($each,' ')+1));
		if(trim(substr($each, -2)) == $plus_o){
			file_put_contents($outputFile.$tail, $contigSeq[$partEach], FILE_APPEND);
		}
		else{
			file_put_contents($outputFile.$tail, $contigSeq["$partEach.reverse"], FILE_APPEND);
		}
		$edge = 0;
	}
}

function UnalignFile($oriFile, $Markers_ori, $scaffoldNum, $outputFile, $refPrefix){
	global $plus_o;
	global $contigSeq;
	$unAlign = '';
	$unAlign_seq = '';
	$Contig_dra = array();
	$parse = fopen($oriFile, 'r');
	while(($seq = fgets($parse)) !== false){
		if($seq[0] == '>'){
			array_push($Contig_dra, $seq);
		}
	}
	fclose($parse);
	$n = 0;
	foreach($Contig_dra as $eachInDraft){
		$flag = $eachInDraft;
		$contigNameInD = strtok(substr($flag, 1), " \n");
		foreach($Markers_ori as $eachInContig){
			$contigNameInC = strtok($eachInContig, " \n");
			if($contigNameInD == $contigNameInC){
				$unmapped = 0;
				break;
			}
			else{
				$unmapped = 1;
			}
		}
		if($unmapped == 1){
			$scaffoldNum++;	
			$unAlign .= ">".$refPrefix."Scaffold_".$scaffoldNum."\n";
			$unAlign .= "$contigNameInD $plus_o\n\n";
			$unAlign_seq .= ">".$refPrefix."Scaffold_".$scaffoldNum."\n";
			$unAlign_seq .= $contigSeq[$contigNameInD]."\n";
		}
	}
	file_put_contents($outputFile, $unAlign, FILE_APPEND);
	file_put_contents("$outputFile.fasta", $unAlign_seq, FILE_APPEND);
}

function path_parse($coords){
	$arr = file($coords, FILE_IGNORE_NEW_LINES);
	for($i = 0; $i < 5; $i++){
		array_shift($arr);
	}
	$element = preg_split('/ +/', trim($arr[0]));
	$num = count($element);
	$lines = array();
	foreach($arr as $k => $line){
		$element = preg_split('/ +/', trim($line));
		$lines[$k][7] = $element[7]; 
		$lines[$k][$num-3] = $element[$num-3]; 
		$lines[$k][$num-2] = $element[$num-2]; 
	}
	$line_num = count($lines);
	$pos_len = 0;
	$neg_len = 0;
	for($i = 0; $i < $line_num; $i++){
		if($lines[$i][$num-3] * $lines[$i][$num-2] > 0){
			$pos_len += $lines[$i][7];
		}
		else{
			$neg_len += $lines[$i][7];
		}
	}
	if($pos_len > $neg_len){
		return 1;
	}
	else{
		return -1;
	}
}

function contigGen($contigSet){
	if($contigSet[0][0] == '>'){
		$contigSet_name = $contigSet[0];
		array_shift($contigSet);
	}
	$contigSet = array_reverse($contigSet);
	foreach($contigSet as $k => $contig){
		list($contig_name, $ori) = preg_split('/ +/', $contig);
		$contigSet[$k] = $contig_name.' '.(1-intval($ori));
	}
	if(isset($contigSet_name))
		array_unshift($contigSet, $contigSet_name);
	return $contigSet;
}

function cmp($a, $b){
	if($a[0] == $b[0]){
		return ($a[1] < $b[1]) ? -1 : 1;
	}
	return ($a[0] < $b[0]) ? -1 : 1;
}

function arr_Size($a, $b){
	$a_size = count($a);
	$b_size = count($b);
	if($a_size == $b_size){
		return 0;
	}
	return ($a_size < $b_size) ? 1 : -1;
}

function checkFASTA($file, $subject){
	global $meg, $meg_h;
	$handle = fopen($file, 'r');
	$contents = fread($handle, 1);
	if($contents != '>' && $contents != ';'){
		fclose($handle);
		$subject = $subject == 't' ? 'Target' : 'Reference';
		pErr("$subject is not a FASTA or multi-FASTA file");
	}
	fclose($handle);
}

function pErr($str){
	exit("ERROR: $str\n");
}

function p($str){
	echo "$str\n";
}

function getMicrotime(){
	list($usec, $sec) = explode(' ', microtime());
	return ((float)$usec + (float)$sec);
}

function checkFormat($file,$format){
	$solve = fopen($file, 'r');
	$contents = fread($solve, 1);
	if($format == '1'){
		if($contents != '>' && $contents != ';'){
			fclose($solve);
			$format='FASTA/FA/fasta/fa';
			pErr("Input $file is not a $format file");
		}
	}
	else{
		if($contents != '@' && $contents != ';'){
			fclose($solve);
			$format='FASTQ/FQ/fastq/fq';
			pErr("Input $file is not $format file");
		}
	}
	fclose($solve);
}

function judgeFasta($file){
	return (pathinfo($file)['extension']=="fasta" || pathinfo($file)['extension']=="FASTA" || pathinfo($file)['extension']=="FA" || pathinfo($file)['extension']=="fa");
}

function judgeFastq($file){
	return (pathinfo($file)['extension']=="fastq" || pathinfo($file)['extension']=="FASTQ" || pathinfo($file)['extension']=="FQ" || pathinfo($file)['extension']=="fq");
}

function scoreContig($samFile,$read_length){
	$count=0;
	$maq=0;
	$total=0;
	$rate=round(1/($read_length),9);
	$arr = file($samFile, FILE_IGNORE_NEW_LINES);
	foreach($arr as $k => $line){
		$element = preg_split('/\s+/is',$line);
		if($element[0]!="@SQ" || $element[0]!="@HD"){
			if($element[6]=="=" && $element[5] !="*" && $element[4]!="*"){
				$score[0]=$score[0]+splitStr($element[5],$read_length)[0];
				$score[1]=$score[1]+splitStr($element[5],$read_length)[1];
				$count=$count+1;
				$maq=$maq+$element[4];
				$total=$total+1;
			}
			else{
				if($rate==0.01){
					$score[0]=$score[0]-1;
				}
				else{
					$score[0]=$score[0]-$rate;
				}
				$unalign=$unalign+1;
				$total=$total+1;
			}
		}
	}
	$maqMean=round($maq/($count*42),9);
	if($read_length==100){
		$scoreMean=round(($score[0]-$score[1]*30)/$total,9);
	}
	elseif($read_length==101){
		$scoreMean=round(($score[0]-$score[1])/$total,9);
	}
	else{
		$scoreMean=round(($score[0]-$score[1]*100)/$total,9);
	}
	return $scoreMean;
}

function splitStr($str,$read_length){
	$count_err=0;
	$arr = preg_split("/([a-zA-Z]+)/", $str, 0, PREG_SPLIT_NO_EMPTY | PREG_SPLIT_DELIM_CAPTURE);
	$len = sizeof($arr);
	$rate=round(1/($read_length),9);
	for($i=0;$i<$len;$i++){
		if($arr[$i]=="M"){
			$score[0]=$score[0]+$rate*intval($arr[$i-1]);
		}
		elseif($arr[$i]=="I"||$arr[$i]=="D"||$arr[$i]=="N"||$arr[$i]=="P"||$arr[$i]=="S"||$arr[$i]=="H"){
			//$score[0]=$score[0]-5*intval($arr[$i-1]);
			if($read_length==100){
				$score[1]=$score[1]+1;
			}
			else{
				$score[1]=$score[1]+intval($arr[$i-1]);
			}
		}
		else continue;
	}
	return $score;
}

function alignment($pair1,$pair2,$contigSet,$outPath,$read_length){
	$outPath1=$outPath.'/mapping';
	if(!@mkdir($outPath1, 0777, true) && !is_dir($outPath1)){
		pErr('Can not create the output directory.');
	}
	$CurrentPath = $outPath1.'/'.basename($contigSet);
	system("bowtie2-build $contigSet $outPath1/index");
	system("bowtie2 -x $outPath1/index -1 $pair1 -2 $pair2 -S $CurrentPath.sam");
	//system("rm -f $CurrentPath/index.*");
	if(is_file("$CurrentPath.sam")){
		$score=scoreContig("$CurrentPath.sam",$read_length);
		return $score;
	}
	else{
		pErr("Bowtie running error!");
	}
}

function filterContig($coordsFile,$targetFile){
	print("$coordsFile\n");
	$arr = file($coordsFile, FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES);
	$temp_target= file($targetFile, FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES);
	$temp_target = arrParse($temp_target);
	$i=0;
	$j=0;
	$flag=1;
	$count_filter=0;
	for($i = 0; $i < 5; $i++){
		array_shift($arr);
	}
	foreach($arr as $k => $line){
		$element = preg_split('/ +/', trim($line));
		if($element[9]=="100.00" && $element[11] <= $element[12] && $element[6]==$element[11]){
			$temp_elem=preg_split('/[\s,;]+/',$element[16]);
			$filterEle=">".$temp_elem[0];
			foreach($temp_target as $k => $line2){
				$temp_elem=preg_split('/[\s,;]+/',$line2);
				$temp_judge=$temp_elem[0];
				if($flag==1){
					if($temp_judge==$filterEle){
						$flag=0;
						$count_filter=$count_filter+1;
						$temp_target=array_diff($temp_target,["$line2"]);
					}
					else{
						continue;
					}
				}
				else{
					$count_filter=$count_filter+1;
					$temp_target=array_diff($temp_target,["$line2"]);
					$flag=1;  
				}
			}
		}
	}
	return $temp_target;
}

function Mergecontig($assembled_fasta,$temp_target){
	$count=0;
	$mergResult=Array();
	$arr = file($assembled_fasta, FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES);
	foreach($arr as $k => $line){
		if(strstr($line,'>')){
			$temp_str=">".$count;
			Array_push($mergResult,$temp_str);
			$count=$count+1;
		}
		else{
			$temp_str=$line;
			Array_push($mergResult,$temp_str);
		}
	}
	foreach($temp_target as $k => $line1){
		if(strpos($line1,'>')!==false){
			$temp_str=">".$count;
			Array_push($mergResult,$temp_str);
			$count=$count+1;
		}
		else{
			$temp_str=$line1;
			Array_push($mergResult,$temp_str);
		}
	}
	$mergResult=arrParse($mergResult);
	return $mergResult;
}

function arrTofile($ouputFile,$arr){
	if(is_file($ouputFile)){
		unlink($ouputFile);
	}
	foreach($arr as $k => $line3){
		file_put_contents($ouputFile,"$line3\n",FILE_APPEND);
	}
	if(is_file($ouputFile)){
		return true;
	}
	else{
		return false;
	}
}

function arrParse($arry){
	$count=0;
	$input_str="";
	$temp_str="";
	$arry_new=Array();
	$Last=end($arry);
	foreach($arry as $k => $line){
		if(strstr($line,'>')){
				Array_push($arry_new,$input_str);
				Array_push($arry_new,">".$count);
				$temp_str="";
				$input_str="";
				$count=$count+1;
		}
		else{
			$temp_str=str_replace("\n","",$line);
			$input_str=$input_str.$temp_str;
			if($Last===$line){
				Array_push($arry_new,$input_str);
			}
		}
	}
	$arry_new=array_filter($arry_new);
	return $arry_new;
}

function fileParse($out,$file){
	$arr = file($file, FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES);
	$arr = arrParse($arr);
	if (arrTofile($out,$arr)){
		return $out;
	}
	else{
		pErr("File parsing error!");
	}
}

function my_sort($arrays,$sort_key,$sort_order=SORT_ASC,$sort_type=SORT_NUMERIC){  
	if(is_array($arrays)){
		foreach ($arrays as $array){
			if(is_array($array)){
				$key_arrays[] = $array[$sort_key];
			}
			else{
				return false;
			}
		}
	}
	else{
		return false;
	}
	array_multisort($key_arrays,$sort_order,$sort_type,$arrays);
	return $arrays;
	}

function trimSeq($seqFile1,$seqFile2){
	$n=0;
	$arr=file($seqFile1, FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES);
	$arr_result=Array();
	foreach($arr as $k => $line){
		$l=str_replace("\n","",$line);
		$n=$n+1;
		if($n%4==1){
			Array_push($arr_result,$l);
		}
		elseif($n%4==2){
			if(strlen($l)>=100){
				Array_push($arr_result,substr($l,0,100));
			}
			else{
				Array_push($arr_result,$l);
			}
		}
		elseif($n%4==3){
			Array_push($arr_result,$l);
		}
		else{
			if(strlen($l)>=100){
				Array_push($arr_result,substr($l,0,100));
			}
			else{
				Array_push($arr_result,$l);
			}
		}
	}
	if(arrTofile($seqFile2,$arr_result)){
		return true;
	}
	else{
		pErr('Parsing fastq files error!');
	}
}

////////////Start processing!/////////////////////////////
ini_set('memory_limit', -1);
$totalTimeStart = getMicrotime();
$opt = getopt("e:h", array('eva'));
$meg = 'Usage: php MAC.php <contig_1.fasta> <contig_2.fasta> ... <paired_reads1.fastq> <paired_reads2.fastq> <read length> <insert size> ... <output>';
$meg_h = 'Using option \'-h\' for help';

if(isset($opt['h'])){
	print <<<END
Option:
	<contigs.fasta>         Contig files generated by different assembly methods. At least two contig files are needed, please separate each file with space

	<paired_reads1.fastq>   The left end of Paired-end/mate-pair reads file

	<paired_reads2.fastq>   The right end of Paired-end/mate-pair reads file

	<read length>   The average length of paired reads

	<insert size>   The insert size of paired reads

	<output>   The path to store the result files

	-h         Show help message 

	-e 		Use evaluation

END;
	exit();
}
print("number of arg: $argc\n");
for ($j=0; $j<$argc; $j++){
	if(judgeFasta($argv[$j])){
		$contig_num=$contig_num+1;
	}
}
if($contig_num<2){
	pErr("At least two contig files are needed!");
}
if($argv[$argc-1]=="-e"){
	print("Start evaluation!\n");
	if(($argc-$contig_num-3)%3==0){
		$pair_num=($argc-$contig_num-3)/3;
	}
	else{
		pErr("The Paired-end/mate-pair reads file should be input in pairs and the read length is aldo needed!");
	}

print("Input $contig_num contig files, $pair_num pair files, and start processing\n");
print("====================================================\n");
$CWD = getcwd();
$outPath = $argv[$argc-2] ? $argv[$argc-2] : $CWD.'/Mix_out';
/*
for($m=0; $m< $pair_num; $m++){                                //Obtain fastq file from input
	$fastq_temp="$fastqSet".$m;
	$$fastq_temp=$argv[$contig_num+$m+1];
	if($$fastq_temp != '') checkFormat($$fastq_temp,"2");
}*/
for($f=0;$f<$contig_num;$f++){                                            //Score the contigs set
	$score=0;
	for($m=0; $m< $pair_num; $m++){
		$fastq_temp_1=$contig_num+$m*3+1;
		$fastq_temp_2=$contig_num+$m*3+2;
		$length=$contig_num+$m*3+3;
		$score=alignment($argv[$fastq_temp_1],$argv[$fastq_temp_2],$argv[$f+1],$outPath,$argv[$length])+$score;
		print("score: $score\n");
	}
	$scoreContig[$f]["num"]=$f;
	$scoreContig[$f]["val"]=round($score/$pair_num,7);
	$contig_name=basename($argv[$f+1]);
	print("Aligning $contig_name finished!\n");
}

print_r($scoreContig);
$contigSort=my_sort($scoreContig,"val",SORT_ASC,SORT_NUMERIC);            //Sort the contigs set

print_r($contigSort);
print("Evaluation finished!\n");

for($f=0;$f<$contig_num;$f++){
	$sort=$contigSort[$f]["num"]+1;
	print("$argv[$sort]\n");
}

}
else{
	print("no evaluation\n");
	$outPath = $argv[$argc-1] ? $argv[$argc-1] : $CWD.'/Mix_out';
	$length=$argv[$contig_num+$m*3+3];
	$pair1=$argv[$contig_num+1];
	$pair2=$argv[$contig_num+2];
	$inserSize=$argv[$argc-2];
}


////////////////////Into the main rotate////////////////////////////
for($b=0;$b<$contig_num-1;$b++){
if($argv[$argc-1]=="-e"){
	if($b==0){
		$temp_sort1=$contigSort[0]["num"]+1;
		$temp_sort2=$contigSort[1]["num"]+1;
		$contigSet_1 = $argv[$temp_sort1];
		$contigSet_2 = $argv[$temp_sort2];
	}
	else{
		$contigSet_2= "$outPath/Mixcontig".($b-1).".fasta";
		$temp=$b+1;
		$temp_sort=$contigSort[$temp]["num"];
		$contigSet_1 = $argv[$temp_sort+1];
	}
}
else{
	if($b==0){
		$contigSet_1=$argv[1];
		$contigSet_2=$argv[2];
	}
	else{
		$temp=$b+2;
		$contigSet_1=$argv[$temp];
		$contigSet_2="$outPath/Mixcontig".($b-1).".fasta";
	}
	
}
print("contig_1:$contigSet_1\n");
print("contig_2:$contigSet_2\n");
if($contigSet1 != '') checkFormat($contigSet1,"1");
if($contigSet2 != '') checkFormat($contigSet2,"1");

$mum = 'nuc';
$getRunTime =0;

if(!@mkdir($outPath, 0777, true) && !is_dir($outPath)){
	pErr('Can not create the output directory.');
}

elseif(isset($useCoords)){
	$mummerRunTime = 0;
	$mumOption = "";
}
elseif($mum != ''){
	$refCurrentPath = $outPath.'/'.basename($contigSet_2);
	$mumOption = '';
	$mumOptSuf = $mum;
	$deltaPara = '-r -q'; 
	$mummerRunTime = 0;
	$mummerStartTime = getMicrotime();
	system($mum."mer $mumOption $contigSet_2 $contigSet_1 -p $refCurrentPath.".$mumOptSuf);
	system("delta-filter $deltaPara $refCurrentPath.$mumOptSuf.delta > $refCurrentPath.$mumOptSuf.filter"); 
	system("show-coords -l -d $refCurrentPath.$mumOptSuf.filter > $refCurrentPath.$mumOptSuf.coords");
	$coords = "$refCurrentPath.$mumOptSuf.coords";
	$mummerRunTime = round(getMicrotime() - $mummerStartTime, 3);
}
else{
	pErr($meg_h);
}

if(!file_exists($coords) || $coords == ''){
	pErr($meg_h);
}
//================================================//
$arr = file($coords, FILE_IGNORE_NEW_LINES);                                    //Parsing coords file
$mum = ($arr[1] == 'PROMER') ? 'pro' : 'nuc';
$mumOptSuf = $mum;
for($i = 0; $i < 5; $i++){
	array_shift($arr);
}
$lines = array();
$element = preg_split('/ +/', trim($arr[0]));
$num = count($element);
$weight = 0;
$getMarkerTime = array();
$getMarkerTime[0] = getMicrotime();
$contig1_tag="";
$contig2_tag="";
foreach($arr as $k => $line){
	$element = preg_split('/ +/', trim($line));
	list($contig1_tag, $contig2_tag) = preg_split("/\t/", $element[$num-1]);
	$lines[$contig1_tag][$k][0] = $element[0]; 
	$lines[$contig1_tag][$k][1] = $element[1]; 
	if($lines[$contig1_tag][$k][0] > $lines[$contig1_tag][$k][1]){
		exchange($lines[$contig1_tag][$k][0], $lines[$contig1_tag][$k][1]);	
	}
	$lines[$contig1_tag][$k][3] = $element[3]; 
	$lines[$contig1_tag][$k][4] = $element[4]; 
	if($lines[$contig1_tag][$k][3] > $lines[$contig1_tag][$k][4]){
		exchange($lines[$contig1_tag][$k][3], $lines[$contig1_tag][$k][4]);	
	}
	$lines[$contig1_tag][$k][6] = $element[6]; 
	$lines[$contig1_tag][$k][7] = $element[7];
	$lines[$contig1_tag][$k][$num-3] = $element[$num-3]; 
	$lines[$contig1_tag][$k][$num-2] = $element[$num-2]; 
	$lines[$contig1_tag][$k][$num-1] = $contig2_tag;

}
$getMarkerTime[0] = round(getMicrotime() - $getMarkerTime[0], 3);
$getMarkerTime[1] = getMicrotime();

foreach($lines as $contig_tag => $each_contig){
	usort($lines[$contig_tag], 'cmp');
}

$getMarkerTime[1] = round(getMicrotime() - $getMarkerTime[1], 3);
$getMarkerTime[2] = getMicrotime();
$n = 0;
$ref_marker = '';
foreach($lines as $contig_tag => $each_contig){
	$ref_marker .= ">$contig_tag\n";	
	foreach($each_contig as $i => $each_marker1){
		$n++;
		$ref_marker .= $n."\n";	
		$lines[$contig_tag][$i]['marker'] = $n;
	}
	$ref_marker .= "\n";	
}
$getMarkerTime[2] = round(getMicrotime() - $getMarkerTime[2], 3);
$getMarkerTime[3] = getMicrotime();
$lines2 = array();
$n = 0;
foreach($lines as $contig_tag => $each_contig){
	foreach($each_contig as $i => $each_marker1){
		$contig2_tag = $lines[$contig_tag][$i][$num-1];
		$lines2[$contig2_tag][$n][0] = $lines[$contig_tag][$i][3]; 
		$lines2[$contig2_tag][$n][1] = $lines[$contig_tag][$i][4];
		if($lines[$contig_tag][$i][$num-3] * $lines[$contig_tag][$i][$num-2] > 0){
			$lines2[$contig2_tag][$n]['marker'] = $lines[$contig_tag][$i]['marker'];
		}
		else{
			$lines2[$contig2_tag][$n]['marker'] = -$lines[$contig_tag][$i]['marker']; 
		}
		$n++;
	}
}
$getMarkerTime[3] = round(getMicrotime() - $getMarkerTime[3], 3);
$getMarkerTime[4] = getMicrotime();

foreach($lines2 as $contig_tag => $each_contig){
	usort($lines2[$contig_tag], 'cmp');
}
$getMarkerTime[4] = round(getMicrotime() - $getMarkerTime[4], 3);

$getMarkerTime[5] = getMicrotime();

$target_marker = '';
foreach($lines2 as $contig_tag => $each_contig){
	$target_marker .= ">$contig_tag\n";	
	foreach($each_contig as $i => $each_marker1){
		$target_marker .= $lines2[$contig_tag][$i]['marker']."\n";
	}
	$target_marker .= "\n";	
}
$getMarkerTime[5] = round(getMicrotime() - $getMarkerTime[5], 3);

$coords = basename($coords);

//======================================//
$contigSet = array();

$contigSet[0] = explode('>', $ref_marker);
$contigSet[1] = explode('>', $target_marker);
$ref_originMarkers = $contigSet[0];
$tar_originMarkers = $contigSet[1];

foreach($contigSet as $k => $set){
	$contigNum[$k] = 0;
	$contigs[$k] = array();
	$contigTag[$k] = array();
	$telomere[$k] = array();
	$telomere_label[$k] = array();

	$contig_idx = 0;
	foreach($set as $l => $contig){
		if($contig == '') continue;
		$contigNum[$k]++;
		$markers = explode("\n", trim($contig));
		$markersNum = count($markers);
		$contigTag[$k][$markers[1]] = $markers[0];
		$contigTag[$k][$markers[$markersNum-1]] = $markers[0];
		array_push($telomere[$k], $markers[1]);
		array_push($telomere[$k], $markers[$markersNum-1] * -1);
		$telomere_label[$k][$markers[1]] = $contig_idx;
		$telomere_label[$k][$markers[$markersNum-1] * -1] = $contig_idx;
		$contig_idx++;
		array_shift($markers);
		array_push($contigs[$k], $markers);
	}
}

$contigSet = $contigs;

$stepTime = array();
$MAC_StartTime = getMicrotime();
//======================================//
$stepTime[0] = getMicrotime();
foreach($contigSet as $k => $set){
	foreach($set as $l => $contig){
		$markersNum = count($contig);
		for($i = 0; $i < $markersNum; $i++){
			array_push($contigSet[$k][$l], $contig[$markersNum-1 - $i] * -1); 
		}
	}
	
}
$stepTime[0] = round(getMicrotime() - $stepTime[0], 3);

//======================================//
$stepTime[1] = getMicrotime();
unset($reverse_contig);
foreach($contigSet[0] as $l => $contig){
	$reverse_contig[$l] = array_reverse($contig);
}

print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
$sigma_piInverse = multiply($contigSet[1], $reverse_contig);

$index_in_sigma_piI = array();
foreach($sigma_piInverse as $key => $cycle){
	if(count($cycle) < 2) continue;
	foreach($cycle as $marker){
		$index_in_sigma_piI[$marker] = $key;
	}
}

$stepTime[1] = round(getMicrotime() - $stepTime[1], 3);

//======================================//
$stepTime[2] = getMicrotime();
if(0)
foreach($telomere[0] as $i => $x){
	if(!isset($index_in_sigma_piI[$x])) continue; 
	for($j = $i%2==0? $i+2: $i+1; $j<$contigNum[0]*2; $j++){
		if(!isset($telomere[0][$j])) continue;
		$y = $telomere[0][$j];
		if(!isset($index_in_sigma_piI[$y])) continue; 
		if($index_in_sigma_piI[$x] == $index_in_sigma_piI[$y] && $telomere_label[0][$x] != $telomere_label[0][$y]){
			$contigSet[0] = connect(array($x, $y), $contigSet[0], $telomere_label[0][$x], $telomere_label[0][$y]);
			unset($telomere[0][$i]); 
			unset($telomere[0][$j]); 
			vertex_parse($telomere_label[0], $x, $y);
		}
	}
}

$contigNum_temp = $contigNum[1] << 1; 
foreach($telomere[1] as $i => $x){
	if(!isset($index_in_sigma_piI[$x])) continue; 
	for($j = $i%2==0? $i+2: $i+1; $j<$contigNum_temp; $j++){ 
		if(!isset($telomere[1][$j])) continue;
		$y = $telomere[1][$j];
		if(!isset($index_in_sigma_piI[$y])) continue; 
		if($index_in_sigma_piI[$x] == $index_in_sigma_piI[$y] && $telomere_label[1][$x] != $telomere_label[1][$y]){
			$contigSet[1] = connect(array($x, $y), $contigSet[1], $telomere_label[1][$x], $telomere_label[1][$y]);
			unset($telomere[1][$i]); 
			unset($telomere[1][$j]); 
			vertex_parse($telomere_label[1], $x, $y);
		}
	}
}

$contigSet[1] = array_values($contigSet[1]);
$stepTime[2] = round(getMicrotime() - $stepTime[2], 3);
$MAC_RunTime = round(getMicrotime() - $MAC_StartTime, 3);

//======================================//
$stepTime[3] = getMicrotime();
$plus_strand = array();
foreach($contigSet[1] as $j => $eachContig){
	$markerNum = count($eachContig) / 2;
	foreach($eachContig as $k => $marker){
		foreach($telomere[1] as $eachTelomere){
			if($marker == $eachTelomere){
				$plus_strand[1][$j] = cycle($eachContig, $k);
				array_splice($plus_strand[1][$j], $markerNum);
				break(2);
			}
		}
	}
}

$plus_o = '0';
$minus_o = '1';
$scaffoldNum = 0;

$temp = '';
unset($scaffold);
foreach($plus_strand[1] as $j => $eachContig){
	$k = $j+1;
	$major_ori = 0;
	$scaffold[$k] = array();
	foreach($eachContig as $marker){
		@$tagTmp = $contigTag[1][$marker];
		@$tagTmpN = $contigTag[1][$marker*-1];
		if(isset($tagTmp) && $tagTmp != $temp){
			array_push($scaffold[$k], "$tagTmp $plus_o");
			$temp = $tagTmp;
			$major_ori++;
		}
		else if(isset($tagTmpN) && $tagTmpN != $temp){
			array_push($scaffold[$k], "$tagTmpN $minus_o");
			$temp = $tagTmpN;
			$major_ori--;
		}
	}
	if($major_ori < 0){
		$scaffold[$k] = contigGen($scaffold[$k]);
	}
}
$scaffoldNum = $j+1;

usort($scaffold, 'arr_Size');

$all_scaffolds = '';
foreach($scaffold as $i => $scaffold_each){
	$all_scaffolds .= '>Scaffold_'.($i+1)."\n";
	foreach($scaffold_each as $line){
		$all_scaffolds .= "$line\n";
	}
	$all_scaffolds .= "\n";
}

$outputFile_target = "$outPath/Mixcontig".$b;
file_put_contents($outputFile_target, $all_scaffolds);

$stepTime[3] = round(getMicrotime() - $stepTime[3], 3);
$genSeqTime = getMicrotime();

$contigSeq = array();
if($contigSet_1 != ''){
	parse_contig($contigSet_1);
	$fastaFile = 1;
	if($fastaFile == 1){
		fastaFile($contigSet_1, $outputFile_target);
	}
	$genUnmapped = 1;
	if($genUnmapped == 1){
		UnalignFile($contigSet_1, $tar_originMarkers, $scaffoldNum, $outputFile_target, '');
	}
	$genSeqTime = round((getMicrotime() - $genSeqTime), 3);
		$dotplotStartTime = getMicrotime();
		$assembFile = '';
		$draft_assembly = $outputFile_target;
		$assembled_fasta = "$draft_assembly.fasta";

		
}

$totalTime = round((getMicrotime() - $totalTimeStart), 3);
if($getRunTime == 1){
	$runtimeContext = "Runtime info:\n";
	$runtimeContext .= "MUMmer: $mummerRunTime seconds\n";
	$runtimeContext .= "Main algorithm: $MAC_RunTime seconds\n";
	$runtimeContext .= "Generate the fasta file: $genSeqTime seconds\n";
	$runtimeContext .= "Total Time: $totalTime seconds\n";
	p($runtimeContext);
	file_put_contents("$outPath/MAC_$mum.runTime", $runtimeContext);
}

print("Merging finished!\n");
}

if(is_file("$outputFile_target.fasta")){
	$mumOption="";
	$deltaPara='-r -q';
	if($argv[$argc-1]=="-e"){
		$sort=$contigSort[$contig_num-1]["num"]+1;
	}
	else{
		$sort=$contig_num;
	}
	$assembled_fasta = "$outPath/Mixcontig".($b-1).".fasta";
	$parse_arg=fileParse("$outPath/parsefile.fasta",$argv[$sort]);
	print("$parse_arg\n");
	print("$assembled_fasta\n");
	system("nucmer $mumOption $parse_arg $mumOption $assembled_fasta -p $outPath/Mix2");
	system("delta-filter $deltaPara $outPath/Mix2.delta > $outPath/Mix2.filter");
	system("show-coords -l -d $outPath/Mix2.filter > $outPath/Mix2.coords");
	
	$temp_target=filterContig("$outPath/Mix2.coords","$argv[$sort]");
	$mergResult=Mergecontig("$assembled_fasta",$temp_target);

	if(arrTofile("$outPath/MixOut.fasta",$mergResult)){
		print("====================================================\n");
		print("Mix running finished!\n");
	}
	else{
		pErr("Mix running error!\n");
	}
	
}
else{
	pErr("Mix running error!\n");
}


?>
