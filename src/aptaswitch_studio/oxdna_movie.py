"""Offline playback/export of actual coarse-grained oxDNA trajectories."""
from __future__ import annotations

from html import escape
import json
import math
import re
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from aptaswitch_core.localization import tr
from .figure_ui import ExportOption, render_export_menu
from .molecular_viewer_assets import bundled_viewer_javascript


def _contained(root: Path, filename: str) -> Path:
    path = root / filename
    if not filename or path.is_symlink() or not path.resolve().is_relative_to(root.resolve()) or not path.is_file():
        raise ValueError(tr("Fichier de film oxDNA absent ou extérieur au calcul."))
    for parent in path.parents:
        if parent == root:
            break
        if parent.is_symlink():
            raise ValueError(tr("Fichier de film oxDNA absent ou extérieur au calcul."))
    return path


def read_movie(topology: Path, trajectory: Path, *, max_frames=600) -> dict:
    """Validate new-format topology and stream frame offsets before subsampling."""
    lines = topology.read_text().splitlines()
    header = lines[0].split()
    if len(header) != 3 or header[2] != "5->3" or int(header[1]) != 2:
        raise ValueError(tr("Topologie du film oxDNA invalide."))
    sequence, strands, chemistries = [], [], []
    for index, line in enumerate(lines[1:]):
        fields = line.split()
        if not fields:
            continue
        bases = fields[0]
        chemistry = next((field.split("=", 1)[1] for field in fields[1:] if field.startswith("type=")), "DNA")
        if set(bases) - set("ACGTU") or chemistry not in {"DNA", "RNA"}:
            raise ValueError(tr("Topologie du film oxDNA invalide."))
        sequence.extend(bases)
        strands.extend([index] * len(bases))
        chemistries.extend([chemistry] * len(bases))
    n = len(sequence)
    if not 2 <= n <= 5000 or n != int(header[0]) or set(strands) != {0, 1}:
        raise ValueError(tr("Topologie du film oxDNA invalide."))
    if trajectory.stat().st_size > 512 * 1024 * 1024:
        raise ValueError(tr("La trajectoire est trop volumineuse pour le lecteur intégré ; téléchargez-la pour oxView."))
    offsets = []
    with trajectory.open() as stream:
        while True:
            start = stream.tell()
            line = stream.readline()
            if not line:
                break
            if not line.strip():
                continue
            if not line.startswith("t ="):
                raise ValueError(tr("Trame oxDNA invalide."))
            for _ in range(n + 2):
                if not stream.readline():
                    raise ValueError(tr("La dernière trame oxDNA est incomplète."))
            offsets.append(start)
    if not offsets:
        raise ValueError(tr("Aucune trame disponible dans cette trajectoire."))
    count = min(max_frames, len(offsets), max(2, 120000 // n))
    selected = sorted({round(i * (len(offsets) - 1) / max(1, count - 1)) for i in range(count)})
    frames, steps, boxes = [], [], []
    with trajectory.open() as stream:
        for index in selected:
            stream.seek(offsets[index])
            step = float(stream.readline().split("=", 1)[1])
            if not math.isfinite(step):
                raise ValueError(tr("Trame oxDNA invalide."))
            steps.append(step)
            box = [float(value) for value in stream.readline().split("=", 1)[1].split()]
            if len(box) != 3 or not all(math.isfinite(v) and v > 0 for v in box):
                raise ValueError(tr("Boîte oxDNA invalide."))
            stream.readline()
            frame = []
            for _ in range(n):
                values = [float(value) for value in stream.readline().split()]
                if len(values) != 15 or not all(math.isfinite(v) for v in values):
                    raise ValueError(tr("Coordonnées oxDNA invalides."))
                frame.append([round(v, 5) for v in values[:9]])
            frames.append(frame)
            boxes.append(box)
    return {"sequence": sequence, "strands": strands, "chemistries": chemistries, "frames": frames, "steps": steps,
            "boxes": boxes, "source_frames": len(offsets), "selected_frames": selected}


def apply_domain_colors(movie: dict, visualization: dict | None) -> dict:
    """Color only a sequence-matched topology, identically in live and film views."""
    from aptaswitch_core.preview import DOMAIN_COLORS, TRIGGER_COLOR
    result = dict(movie)
    if visualization:
        expected = "".join(visualization.get("sequence", [])).upper().replace("U", "T")
        actual = "".join(movie["sequence"]).upper().replace("U", "T")
        colors, domains = visualization.get("colors", []), visualization.get("domains", [])
        if (expected != actual or len(colors) != len(actual) or len(domains) != len(actual)
                or any(not isinstance(color, str) or not re.fullmatch(r"#[0-9a-fA-F]{6}", color) for color in colors)
                or any(not isinstance(label, str) for label in domains)):
            raise ValueError(tr("Les couleurs des régions ne correspondent pas aux séquences du film."))
        result.update(colors=list(colors), domains=list(domains))
    else:
        result["colors"] = [TRIGGER_COLOR if strand else DOMAIN_COLORS["Stem L"] for strand in movie["strands"]]
        result["domains"] = ["Trigger" if strand else "Switch" for strand in movie["strands"]]
    result["legend"] = [(tr(label), color) for label, color in dict.fromkeys(zip(result["domains"], result["colors"]))]
    return result


def _legend_groups(movie: dict) -> list[tuple[str, list[tuple[str, str]]]]:
    """Use the topology to keep trigger extensions with their own strand."""
    groups = []
    for strand, label in ((0, movie.get("first_strand_label", "Toehold switch")), (1, "Trigger")):
        entries = dict.fromkeys((tr(domain), color) for index, domain, color in
                               zip(movie["strands"], movie["domains"], movie["colors"]) if index == strand)
        if entries:
            groups.append((tr(label), list(entries)))
    return groups


def _legend_html(movie: dict) -> str:
    groups = []
    for label, entries in _legend_groups(movie):
        chips = "".join(
            f'<span class="oxdna-domain"><i style="background:{escape(color, quote=True)}" aria-hidden="true"></i>{escape(domain)}</span>'
            for domain, color in entries
        )
        groups.append(f'<div class="oxdna-legend-group"><strong>{escape(label)}</strong><div class="oxdna-domains">{chips}</div></div>')
    return """<style>
.oxdna-legend{font-family:system-ui,sans-serif;color:#233147;padding:14px 16px;border:1px solid #e2e8f0;border-radius:12px;background:#f8fafc}
.oxdna-legend-title{font-size:13px;font-weight:600;margin-bottom:12px}
.oxdna-legend-group{display:grid;grid-template-columns:120px minmax(0,1fr);gap:10px;margin-top:10px;font-size:13px}
.oxdna-legend-group strong{font-size:12px;padding-top:4px}.oxdna-domains{display:flex;flex-wrap:wrap;gap:6px 12px}
.oxdna-domain{display:inline-flex;align-items:center;gap:6px;line-height:1.5}.oxdna-domain i{width:10px;height:10px;border-radius:3px;flex-shrink:0;border:1px solid #0001}
@media(max-width:520px){.oxdna-legend-group{grid-template-columns:1fr;gap:4px}}
</style>""" + f'<section class="oxdna-legend" aria-label="{escape(tr("Légende des régions"), quote=True)}"><div class="oxdna-legend-title">{escape(tr("Légende des régions"))}</div>{"".join(groups)}</section>'


def _movie_title(*, endpoint: str, reached: bool, preparation: bool, live: bool,
                 replica_number: int | None = None) -> str:
    if endpoint == "free_md":
        title = tr("MD libre ADN–ADN · aptamère et trigger")
    elif endpoint == "first_intended_base_pair":
        title = tr("Trajectoire FFS · premier appariement cible atteint") if reached else tr("Trajectoire FFS partielle · premier appariement cible non atteint")
    elif endpoint == "opening":
        title = tr("Trajectoire FFS ayant atteint le critère d’ouverture") if reached else tr("Trajectoire FFS partielle — ouverture non atteinte")
    else:
        title = tr("Trajectoire FFS · critère final non identifié")
    if preparation:
        title = tr("Switch et trigger · configuration gros grain préparée")
    if live:
        title = tr("Vue live · MD libre ADN–ADN") if endpoint == "free_md" else tr("Vue 3D en direct")
    if type(replica_number) is int and replica_number > 0:
        title += " · " + tr("Réplique {number}", number=replica_number)
    return title


def movie_html(movie: dict, *, endpoint: str = "opening", reached: bool = False,
               preparation: bool = False, live: bool = False, view_key: str = "",
               replica_number: int | None = None, embedded: bool = False) -> str:
    """Share rendering and controls; standalone exports retain title and legend."""
    if "colors" not in movie:
        movie = apply_domain_colors(movie, None)
    data = json.dumps({**movie, "live": live, "view_key": view_key}, allow_nan=False, separators=(",", ":")).replace("<", "\\u003c")
    title = _movie_title(endpoint=endpoint, reached=reached, preparation=preparation,
                         live=live, replica_number=replica_number)
    labels = json.dumps({"play": tr("Lecture"), "pause": tr("Pause"), "record": tr("Exporter la vidéo"),
                         "recording": tr("Enregistrement…"), "unavailable": tr("Export vidéo indisponible dans ce navigateur ; utilisez le film HTML."),
                         "frame": tr("Trame"), "step": tr("Pas MD")}, ensure_ascii=False).replace("<", "\\u003c")
    return """<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>__TITLE__</title>
<style>
*{box-sizing:border-box}body{font-family:system-ui,sans-serif;margin:0;color:#233147;background:white}h3{font-size:16px;font-weight:600;margin:16px}
main{display:flex;flex-direction:column;min-width:0}#view{width:100%;height:480px;position:relative;min-height:200px;overflow:hidden;border-radius:10px}
#bar{display:flex;align-items:center;gap:8px;padding:10px 4px;flex-wrap:wrap;border-top:1px solid #edf1f5}#playback{display:contents}#playback[hidden]{display:none}
button,select{font:inherit;font-size:13px;border:1px solid #d8e0ea;border-radius:8px;background:#f8fafc;color:#233147;padding:8px 12px;cursor:pointer}
button:hover{background:#edf3fa}button:focus-visible,select:focus-visible,input:focus-visible{outline:2px solid #1258dc;outline-offset:2px}
button:disabled{opacity:.55;cursor:default}input{flex:1;min-width:90px;accent-color:#1258dc}#info{font-size:12px;color:#64748b;padding:0 4px 8px}
.embedded{height:100vh;overflow:hidden}.embedded main{height:100%}.embedded #view{height:auto;flex:1;min-height:0}.embedded #bar,.embedded #info{flex-shrink:0}
.standalone .oxdna-legend{margin:8px 16px 16px}
</style></head><body class="__MODE__"><main>
__HEADER__<div id="view"></div><div id="bar" role="group" aria-label="__CONTROLS__"><button id="fit">__FIT__</button><div id="playback"__HIDE_PLAYBACK__><button id="play"></button><input id="slider" aria-label="__POSITION__" type="range" min="0" value="0"><select id="speed" aria-label="__SPEED__"><option value="6">6 fps</option><option value="12" selected>12 fps</option><option value="24">24 fps</option></select><button id="record"></button></div></div><div id="info"__HIDE_INFO__></div>
__LEGEND__</main>
<script>__VENDOR__</script><script>
const film=__FILM__, L=__LABELS__;
const viewer=$3Dmol.createViewer(document.getElementById('view'),{backgroundColor:'white',antialias:true});
const slider=document.getElementById('slider'), play=document.getElementById('play'), record=document.getElementById('record');
slider.max=film.frames.length-1; play.textContent=L.play;record.textContent=L.record;
let current=0,timer=null,recording=false,saveCamera=()=>{};
const low=[Infinity,Infinity,Infinity],high=[-Infinity,-Infinity,-Infinity];
const display=film.frames.map((f,index)=>{
 const box=film.boxes[index],centres=[[0,0,0],[0,0,0]],counts=[0,0],positions=[];
 f.forEach((row,i)=>{
  const s=film.strands[i],p=row.slice(0,3);
  if(i>0&&film.strands[i-1]===s)for(let a=0;a<3;a++){const d=row[a]-f[i-1][a];p[a]=positions[i-1][a]+d-box[a]*Math.round(d/box[a]);}
  positions.push(p);counts[s]++;for(let a=0;a<3;a++)centres[s][a]+=p[a];
 });
 centres.forEach((c,s)=>c.forEach((v,a)=>c[a]=v/counts[s]));
 const shift=centres[1].map((v,a)=>-box[a]*Math.round((v-centres[0][a])/box[a]));
 return positions.map((p,i)=>p.map((v,a)=>{const coordinate=(v+(film.strands[i]?shift[a]:0)-centres[0][a])*8.518;low[a]=Math.min(low[a],coordinate-5);high[a]=Math.max(high[a],coordinate+5);return coordinate;}));
});
function fit(){
 const corners=[];for(let i=0;i<8;i++)corners.push({elem:'C',x:(i&1)?high[0]:low[0],y:(i&2)?high[1]:low[1],z:(i&4)?high[2]:low[2]});
 const extent=viewer.addModel();extent.addAtoms(corners);viewer.zoomTo();viewer.removeModel(extent);viewer.__fitted=true;
}
function draw(index){
 current=index;slider.value=index;
 const f=film.frames[index];
 const atoms=[];
 f.forEach((row,i)=>{
  const s=film.strands[i], r=display[index][i],a1=row.slice(3,6),a3=row.slice(6,9);
  const a2=[a3[1]*a1[2]-a3[2]*a1[1],a3[2]*a1[0]-a3[0]*a1[2],a3[0]*a1[1]-a3[1]*a1[0]];
  // The actual NA interaction sites: DNA2's groove offset and RNA2's tilt.
  const backbone=r.map((v,a)=>v+8.518*(film.chemistries[i]==='RNA'?-.4*a1[a]+.2*a3[a]:-.34*a1[a]+.3408*a2[a])),base=r.map((v,a)=>v+.4*8.518*a1[a]);
  const bonds=[2*i+1];if(i>0&&film.strands[i-1]===s)bonds.push(2*(i-1));if(i+1<f.length&&film.strands[i+1]===s)bonds.push(2*(i+1));
  atoms.push({elem:'P',atom:'P',x:backbone[0],y:backbone[1],z:backbone[2],chain:s?'B':'A',resi:i+1,bonds,bondOrder:bonds.map(()=>1)});
  atoms.push({elem:'N',atom:'N',x:base[0],y:base[1],z:base[2],chain:s?'B':'A',resi:i+1,bonds:[2*i],bondOrder:[1]});
 });
 viewer.removeAllModels();viewer.addModel().addAtoms(atoms);
 for(const color of new Set(film.colors)){
  const resi=film.colors.flatMap((value,i)=>value===color?[i+1]:[]);
  viewer.setStyle({resi},{stick:{radius:.7,color},sphere:{radius:1.3,color}});
 }
 if(!viewer.__fitted)fit();
 viewer.render();document.getElementById('info').textContent=(film.live?'':L.frame+' '+(index+1)+' / '+film.frames.length+' · ')+L.step+' '+film.steps[index];
}
document.getElementById('fit').onclick=()=>{fit();viewer.render();saveCamera();};
function stop(){if(timer)clearInterval(timer);timer=null;play.textContent=L.play;}
function start(){stop();play.textContent=L.pause;timer=setInterval(()=>draw((current+1)%film.frames.length),1000/Number(document.getElementById('speed').value));}
play.onclick=()=>timer?stop():start();slider.oninput=()=>{stop();draw(Number(slider.value));};document.getElementById('speed').onchange=()=>{if(timer)start();};
record.onclick=async()=>{
 if(recording)return;stop();
 const canvas=document.querySelector('#view canvas');
 if(!canvas.captureStream||!window.MediaRecorder){alert(L.unavailable);return;}
 const type=['video/webm;codecs=vp9','video/webm','video/mp4'].find(t=>MediaRecorder.isTypeSupported(t));
 if(!type){alert(L.unavailable);return;}
 recording=true;record.disabled=true;play.disabled=true;slider.disabled=true;record.textContent=L.recording;
 const fps=Number(document.getElementById('speed').value),stream=canvas.captureStream(fps),chunks=[];
 try{
  const recorder=new MediaRecorder(stream,{mimeType:type});recorder.ondataavailable=e=>{if(e.data.size)chunks.push(e.data);};
  const finished=new Promise(resolve=>{recorder.onstop=resolve;});recorder.start();
  for(let i=0;i<film.frames.length;i++){draw(i);await new Promise(r=>setTimeout(r,1000/fps));}
  recorder.stop();await finished;
  const url=URL.createObjectURL(new Blob(chunks,{type})),a=document.createElement('a');a.href=url;a.download='__VIDEO_NAME__.'+(type.startsWith('video/mp4')?'mp4':'webm');a.click();setTimeout(()=>URL.revokeObjectURL(url),10000);
 }finally{stream.getTracks().forEach(t=>t.stop());recording=false;record.disabled=false;play.disabled=false;slider.disabled=false;record.textContent=L.record;}
};draw(0);
if(film.view_key){
 const storageKey='oxdna-camera-'+film.view_key;
 try{const saved=JSON.parse(sessionStorage.getItem(storageKey));if(Array.isArray(saved)&&saved.length===8&&saved.every(Number.isFinite)){viewer.setView(saved);viewer.render();}}catch(e){}
 saveCamera=()=>{try{sessionStorage.setItem(storageKey,JSON.stringify(viewer.getView()));}catch(e){}};
 saveCamera();
 ['pointerup','touchend','wheel'].forEach(event=>document.getElementById('view').addEventListener(event,()=>setTimeout(saveCamera,80)));
}
new ResizeObserver(()=>{viewer.resize();viewer.render();}).observe(document.getElementById('view'));
</script></body></html>""".replace("__MODE__", "embedded" if embedded else "standalone").replace(
        "__HEADER__", "" if embedded else f"<h3>{escape(title)}</h3>").replace(
        "__LEGEND__", "" if embedded else _legend_html(movie)).replace(
        "__HIDE_PLAYBACK__", " hidden" if preparation or live else "").replace(
        "__HIDE_INFO__", " hidden" if preparation else "").replace(
        "__CONTROLS__", escape(tr("Commandes de la vue 3D"), quote=True)).replace(
        "__FIT__", escape(tr("Recentrer"))).replace(
        "__POSITION__", escape(tr("Position dans le film"), quote=True)).replace(
        "__SPEED__", escape(tr("Vitesse de lecture"), quote=True)).replace(
        "__VIDEO_NAME__", "oxDNA_free_MD" if endpoint == "free_md" else "oxDNA_FFS").replace(
        "__TITLE__", escape(title)).replace("__VENDOR__", bundled_viewer_javascript()).replace("__FILM__", data).replace("__LABELS__", labels)


def render_oxdna_movie(result: dict, run_dir: Path, key: str) -> None:
    filename = result.get("movie_file") or result.get("trajectory_file")
    if not filename:
        return
    try:
        trajectory = _contained(run_dir, filename)
        topology = _contained(run_dir, result.get("topology_file", "system.top"))
        movie = apply_domain_colors(read_movie(topology, trajectory), result.get("visualization"))
        preparation = result.get("movie_kind") == "prepared" or result.get("status") == "prepared"
        live = result.get("movie_kind") == "live"
        endpoint = result.get("trajectory_endpoint") or (result.get("kinetics") or {}).get("endpoint", "opening")
        if endpoint == "free_md":
            movie["first_strand_label"] = "Aptamère"
        options = dict(endpoint=endpoint, reached=result.get("final_interface_reached") is True,
                       preparation=preparation, live=live, replica_number=result.get("replica_number"))
        html = movie_html(movie, **options, view_key=key, embedded=True)
    except (OSError, ValueError, KeyError, IndexError) as exc:
        st.info(tr(str(exc)))
        return
    if endpoint == "free_md":
        explanation = tr("Trajectoire continue de MD libre ADN–ADN. Les couleurs identifient les brins et la zone surveillée ; la vitesse de lecture n’est pas une durée expérimentale.")
    elif live:
        explanation = tr("Vue actualisée pendant le calcul, avec les couleurs de la 2D. La FFS teste plusieurs branches : un retour à une configuration précédente est normal.")
    elif preparation:
        explanation = tr("Configuration gros grain colorée par région, comme la 2D. Les contraintes de préparation ont été retirées avant le contrôle final.")
    else:
        explanation = tr("Film gros grain coloré par région, comme la 2D et la vue live. Il suit une branche FFS réellement simulée ; sa vitesse de lecture n’est pas une durée expérimentale.")
    with st.container(border=True):
        st.markdown("**" + _movie_title(**options) + "**", help=explanation)
        components.html(html, height=540 if live else 500 if preparation else 620, scrolling=False)
    st.html(_legend_html(movie))
    if movie["source_frames"] > len(movie["frames"]):
        st.caption(tr("Le lecteur affiche un sous-ensemble des trames ; la trajectoire téléchargée conserve les trames enregistrées."))
    if live:
        return
    standalone = movie_html(movie, **options)
    render_export_menu((
        ExportOption("Télécharger la vue 3D autonome" if preparation else "Télécharger le film HTML autonome",
                     standalone.encode(), "oxDNA_3D.html" if preparation else "oxDNA_free_MD_film.html" if endpoint == "free_md" else "oxDNA_FFS_film.html", "text/html", key + "_html"),
        ExportOption("Télécharger la configuration oxDNA" if preparation else "Télécharger la trajectoire oxDNA",
                     trajectory.read_bytes(), "configuration.dat" if preparation else "trajectory.dat", "text/plain", key + "_traj"),
        ExportOption("Télécharger la topologie oxDNA", topology.read_bytes(), "system.top", "text/plain", key + "_top"),
    ))
