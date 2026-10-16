import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';
const base=process.env.NCOMMS_REVIEWER_ROOT || path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const input=path.join(base,'REPRODUCED/source_data_payload.json');
const target=process.env.NCOMMS_XLSX_OUTPUT || path.join(base,'REPRODUCED/Source_Data.xlsx');
const preview=path.join(base,'REPRODUCED/workbook_previews');await fs.mkdir(preview,{recursive:true});
const payload=JSON.parse(await fs.readFile(input,'utf8'));const wb=Workbook.create();
const col=n=>{let s='';while(n){n--;s=String.fromCharCode(65+n%26)+s;n=Math.floor(n/26);}return s;};
for(const {name,matrix} of payload){
 const sh=wb.worksheets.add(name), nr=matrix.length,nc=matrix[0].length;
 const r=sh.getRangeByIndexes(0,0,nr,nc);r.values=matrix;
 r.format.font={name:'Arial',size:10,color:'#202020'};r.format.rowHeight=32;
 sh.getRangeByIndexes(0,0,1,nc).format={fill:'#174A67',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,rowHeight:46};
 for(let c=0;c<nc;c++){
  const text=matrix.slice(1).some(r=>typeof r[c]==='string'&&r[c].length>22);
  sh.getRangeByIndexes(0,c,nr,1).format.columnWidth=text?44:21;
  sh.getRangeByIndexes(1,c,nr-1,1).format.wrapText=true;
  if(matrix.slice(1).every(r=>typeof r[c]==='number'))sh.getRangeByIndexes(1,c,nr-1,1).format.numberFormat=matrix[0][c]==='year'?'0':matrix.slice(1).some(r=>!Number.isInteger(r[c]))?'0.0000':'#,##0';
 }
 for(let i=1;i<nr;i++)sh.getRangeByIndexes(i,0,1,nc).format.rowHeight=Math.max(32,...matrix[i].map(v=>typeof v==='string'?Math.ceil(v.length/43)*14+4:32));
 sh.freezePanes.freezeRows(1);sh.freezePanes.freezeColumns(2);sh.showGridLines=false;
 sh.tables.add('A1:'+col(nc)+nr,true,'SD_'+name);
}
wb.recalculate();
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:100},summary:'Source Data error scan'});
await fs.writeFile(path.join(preview,'error_scan.ndjson'),errors.ndjson);
for(const {name,matrix} of payload){
 const png=await wb.render({sheetName:name,range:`A1:${col(Math.min(6,matrix[0].length))}${Math.min(6,matrix.length)}`,scale:1,format:'png'});
 await fs.writeFile(path.join(preview,name+'.png'),new Uint8Array(await png.arrayBuffer()));
}
await (await SpreadsheetFile.exportXlsx(wb)).save(target);console.log('Exported',payload.length,'sheets:',target);
