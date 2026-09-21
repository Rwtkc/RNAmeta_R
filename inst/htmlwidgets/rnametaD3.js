HTMLWidgets.widget({
  name: "rnametaD3",
  type: "output",
  factory: function(el) {
    const WIDTH = 1180, HEIGHT = 560;
    const MARGIN = { top: 92, right: 36, bottom: 108, left: 82 };
    const PANEL = "#fffaf3", GRID = "rgba(133,155,122,0.16)", AXIS = "#41503c", TITLE = "#22301f";
    const list = (value) => value == null ? [] : (Array.isArray(value) ? value : (typeof value === "object" ? Object.values(value) : [value]));
    const finite = (values) => list(values).map(Number).filter(Number.isFinite);
    const axisLabel = (value, percent) => percent ? `${Math.round(value)}%` : (value >= 1000 ? `${(value / 1000).toFixed(1)}k` : String(Math.round(value * 10) / 10));

    function tooltip(tooltipNode, container, event, title, rows) {
      tooltipNode.replaceChildren();
      const heading = document.createElement("div");
      heading.className = "rnameta-d3-tooltip__label"; heading.textContent = title;
      tooltipNode.appendChild(heading);
      rows.forEach(([keyText, valueText]) => {
        const row = document.createElement("div"); row.className = "rnameta-d3-tooltip__row";
        const key = document.createElement("span"); key.textContent = keyText;
        const value = document.createElement("strong"); value.textContent = valueText;
        row.append(key, value); tooltipNode.appendChild(row);
      });
      tooltipNode.dataset.visible = "true";
      const bounds = container.getBoundingClientRect(), tip = tooltipNode.getBoundingClientRect();
      const px = event.clientX - bounds.left, py = event.clientY - bounds.top;
      const right = px + 22, left = px - tip.width - 18;
      tooltipNode.style.left = `${right + tip.width > bounds.width ? Math.max(10, left) : right}px`;
      tooltipNode.style.top = `${Math.min(Math.max(10, py - tip.height / 2), bounds.height - tip.height - 10)}px`;
    }
    const hideTooltip = (node) => { node.dataset.visible = "false"; };

    function baseChart(config) {
      el.replaceChildren(); el.className = "rnameta-d3-widget";
      const container = document.createElement("div"); container.className = "rnameta-d3-chart";
      const svgNode = document.createElementNS("http://www.w3.org/2000/svg", "svg");
      svgNode.setAttribute("viewBox", `0 0 ${WIDTH} ${HEIGHT}`); svgNode.setAttribute("preserveAspectRatio", "xMinYMin meet");
      svgNode.setAttribute("role", "img"); svgNode.setAttribute("aria-label", config.title || "RNAmeta chart");
      const tooltipNode = document.createElement("div"); tooltipNode.className = "rnameta-d3-tooltip"; tooltipNode.dataset.visible = "false";
      container.append(svgNode, tooltipNode); el.appendChild(container);
      return { container, svgNode, svg: d3.select(svgNode), tooltipNode };
    }

    function drawAxes(panel, x, y, plotWidth, plotHeight, config, percent) {
      const ticks = finite(config.yTicks);
      panel.selectAll(".grid").data(ticks.length ? ticks : y.ticks(5)).join("line")
        .attr("x1", 0).attr("x2", plotWidth).attr("y1", y).attr("y2", y).attr("stroke", GRID);
      panel.append("g").call(d3.axisLeft(y).tickValues(ticks.length ? ticks : null).tickFormat((v) => axisLabel(v, percent)).tickSize(6).tickPadding(10))
        .call((g) => g.select(".domain").attr("stroke", AXIS).attr("stroke-width", 1.2))
        .call((g) => g.selectAll(".tick line").attr("stroke", AXIS))
        .call((g) => g.selectAll("text").attr("fill", AXIS).style("font-size", "13px").style("font-weight", 600));
      panel.append("g").attr("transform", `translate(0,${plotHeight})`).call(d3.axisBottom(x).tickSize(0).tickPadding(14))
        .call((g) => g.select(".domain").attr("stroke", AXIS).attr("stroke-width", 1.2))
        .call((g) => g.selectAll("text").attr("fill", AXIS).style("font-size", "13px").style("font-weight", 700));
    }

    function legend(svg, entries, y) {
      if (entries.length < 2) return;
      const widths = entries.map((d) => Math.max(96, String(d.name).length * 8.2 + 42));
      const gap = 22, total = d3.sum(widths) + gap * (entries.length - 1);
      let cursor = Math.max(MARGIN.left, MARGIN.left + (WIDTH - MARGIN.left - MARGIN.right - total) / 2);
      entries.forEach((entry, index) => {
        const group = svg.append("g").attr("transform", `translate(${cursor},${y})`);
        group.append("rect").attr("width", 15).attr("height", 15).attr("rx", 4).attr("fill", entry.color);
        group.append("text").attr("x", 23).attr("y", 7.5).attr("dy", ".35em").attr("fill", AXIS).style("font-size", "13px").style("font-weight", 600).text(entry.name);
        cursor += widths[index] + gap;
      });
    }

    function titles(svg, config, plotHeight) {
      svg.append("text").attr("x", MARGIN.left).attr("y", 80).attr("fill", TITLE).style("font-size", "15px").style("font-weight", 800).text(config.title || "");
      svg.append("text").attr("transform", `translate(26,${MARGIN.top + plotHeight / 2}) rotate(-90)`).attr("fill", AXIS).style("font-size", "15px").style("font-weight", 700).text(config.yLabel || "");
    }

    function drawBars(config) {
      const chart = baseChart(config), categories = list(config.categories), series = list(config.series);
      if (!categories.length || !series.length) return;
      const plotWidth = WIDTH - MARGIN.left - MARGIN.right, plotHeight = HEIGHT - MARGIN.top - MARGIN.bottom;
      const x = d3.scaleBand().domain(categories).range([0, plotWidth]).padding(0.18);
      const stacked = config.type === "stacked_bar", domain = finite(config.yDomain);
      const allValues = series.flatMap((entry) => finite(entry.values));
      const yMax = domain[1] || (stacked ? 100 : Math.max(1, d3.max(allValues) || 1));
      const y = d3.scaleLinear().domain([0, yMax]).nice().range([plotHeight, 0]);
      const panel = chart.svg.append("g").attr("transform", `translate(${MARGIN.left},${MARGIN.top})`);
      panel.append("rect").attr("width", plotWidth).attr("height", plotHeight).attr("rx", 18).attr("fill", PANEL);
      drawAxes(panel, x, y, plotWidth, plotHeight, config, stacked);
      const data = [];
      categories.forEach((category, ci) => {
        let cumulative = 0;
        series.forEach((entry) => {
          const value = Number(list(entry.values)[ci] || 0), y0 = cumulative;
          data.push({ category, sample: entry.name, originalName: entry.originalName || entry.name, color: entry.color, value, y0, y1: stacked ? cumulative + value : value });
          if (stacked) cumulative += value;
        });
      });
      const inner = d3.scaleBand().domain(series.map((d) => d.name)).range([0, x.bandwidth()]).padding(0.18);
      const bars = panel.selectAll(".rnameta-bar").data(data).join("rect").attr("class", "rnameta-bar")
        .attr("x", (d) => x(d.category) + (stacked ? x.bandwidth() * .22 : inner(d.sample)))
        .attr("width", stacked ? x.bandwidth() * .56 : inner.bandwidth()).attr("fill", (d) => d.color).attr("stroke", "#fff").attr("stroke-width", 1.2)
        .attr("y", (d) => config.animate === false ? y(d.y1) : y(stacked ? d.y0 : 0))
        .attr("height", (d) => config.animate === false ? Math.max(0, y(d.y0) - y(d.y1)) : 0)
        .on("mouseenter", function(event, d) { d3.select(this).attr("stroke", AXIS).attr("stroke-width", 1.5); tooltip(chart.tooltipNode, chart.container, event, d.category, [[stacked ? "Exon Type" : "Sample", d.sample], [stacked ? "Percentage" : "Frequency", stacked ? `${d.value.toFixed(2)}%` : String(d.value)]]); })
        .on("mousemove", function(event, d) { tooltip(chart.tooltipNode, chart.container, event, d.category, [[stacked ? "Exon Type" : "Sample", d.sample], [stacked ? "Percentage" : "Frequency", stacked ? `${d.value.toFixed(2)}%` : String(d.value)]]); })
        .on("mouseleave", function() { d3.select(this).attr("stroke", "#fff").attr("stroke-width", 1.2); hideTooltip(chart.tooltipNode); });
      if (config.animate !== false && !window.matchMedia("(prefers-reduced-motion: reduce)").matches) bars.transition().duration(760).delay((_, i) => i * 30).attr("y", (d) => y(d.y1)).attr("height", (d) => Math.max(0, y(d.y0) - y(d.y1)));
      else bars.attr("y", (d) => y(d.y1)).attr("height", (d) => Math.max(0, y(d.y0) - y(d.y1)));
      if (!stacked && config.showLabels) panel.selectAll(".bar-label").data(data).join("text").attr("class", "bar-label").attr("x", (d) => x(d.category) + inner(d.sample) + inner.bandwidth() / 2).attr("y", (d) => y(d.value) - 8).attr("text-anchor", "middle").attr("fill", AXIS).style("font-size", "12px").style("font-weight", 700).text((d) => d.value > 0 ? d.value : "");
      panel.append("line").attr("x1", 0).attr("x2", plotWidth).attr("y1", plotHeight).attr("y2", plotHeight).attr("stroke", AXIS).attr("stroke-width", 1.2);
      titles(chart.svg, config, plotHeight); legend(chart.svg, series, HEIGHT - 72);
    }

    function summarize(values, transform) {
      const raw = finite(values).filter((v) => transform !== "log2" || v > 0), data = raw.map((v) => transform === "log2" ? Math.log2(v) : v).sort(d3.ascending);
      if (!data.length) return null;
      const q1 = d3.quantileSorted(data, .25), median = d3.quantileSorted(data, .5), q3 = d3.quantileSorted(data, .75), iqr = q3 - q1;
      const inside = data.filter((v) => v >= q1 - 1.5 * iqr && v <= q3 + 1.5 * iqr);
      return { count: data.length, q1, median, q3, lower: d3.min(inside), upper: d3.max(inside), outliers: data.filter((v) => v < q1 - 1.5 * iqr || v > q3 + 1.5 * iqr) };
    }

    function drawBoxPanel(root, groups, width, height, transform, chart, showY) {
      const entries = list(groups).map((group) => ({ ...group, summary: summarize(group.values, transform) })).filter((d) => d.summary);
      if (!entries.length) return;
      const ext = entries.flatMap((d) => [d.summary.lower, d.summary.upper, ...d.summary.outliers]);
      const y = d3.scaleLinear().domain(d3.extent(ext)).nice().range([height, 0]);
      const x = d3.scaleBand().domain(entries.map((d) => d.name)).range([0, width]).padding(.24);
      root.append("rect").attr("width", width).attr("height", height).attr("rx", 16).attr("fill", PANEL);
      root.selectAll(".grid").data(y.ticks(5)).join("line").attr("x1", 0).attr("x2", width).attr("y1", y).attr("y2", y).attr("stroke", GRID);
      if (showY) root.append("g").call(d3.axisLeft(y).ticks(5).tickFormat((v) => transform === "log2" ? axisLabel(2 ** v, false) : axisLabel(v, false)).tickSize(6).tickPadding(10)).call((g) => g.select(".domain").attr("stroke", AXIS).attr("stroke-width", 1.2)).call((g) => g.selectAll(".tick line").attr("stroke", AXIS).attr("stroke-width", 1.1)).call((g) => g.selectAll("text").attr("fill", AXIS).style("font-size", "13px").style("font-weight", 600));
      root.append("g").attr("transform", `translate(0,${height})`).call(d3.axisBottom(x).tickSize(0).tickPadding(14)).call((g) => g.select(".domain").attr("stroke", AXIS).attr("stroke-width", 1.2)).call((g) => g.selectAll("text").attr("fill", AXIS).style("font-size", "12px").style("font-weight", 700));
      const item = root.selectAll(".box-group").data(entries).join("g").attr("class", "box-group").attr("transform", (d) => `translate(${x(d.name)},0)`);
      item.append("line").attr("x1", x.bandwidth()/2).attr("x2", x.bandwidth()/2).attr("y1", (d) => y(d.summary.lower)).attr("y2", (d) => y(d.summary.upper)).attr("stroke", AXIS).attr("stroke-width", 1.6);
      item.append("rect").attr("x", x.bandwidth()*.19).attr("width", x.bandwidth()*.62).attr("y", (d) => y(d.summary.q3)).attr("height", (d) => Math.max(1, y(d.summary.q1)-y(d.summary.q3))).attr("fill", (d) => d.color).attr("fill-opacity", .86).attr("stroke", "#fff").attr("stroke-width", 1.2);
      item.append("line").attr("x1", x.bandwidth()*.19).attr("x2", x.bandwidth()*.81).attr("y1", (d) => y(d.summary.median)).attr("y2", (d) => y(d.summary.median)).attr("stroke", "#fff").attr("stroke-width", 2);
      item.selectAll(".whisker").data((d) => [d.summary.lower,d.summary.upper].map((v)=>({parent:d,value:v}))).join("line").attr("x1",x.bandwidth()*.31).attr("x2",x.bandwidth()*.69).attr("y1",(d)=>y(d.value)).attr("y2",(d)=>y(d.value)).attr("stroke",AXIS).attr("stroke-width",1.6);
      item.selectAll(".outlier").data((d)=>d.summary.outliers.map((v)=>({parent:d,value:v}))).join("circle").attr("cx",x.bandwidth()/2).attr("cy",(d)=>y(d.value)).attr("r",3.5).attr("fill",(d)=>d.parent.color).attr("stroke","#111");
      item.append("rect").attr("width",x.bandwidth()).attr("height",height).attr("fill","transparent").style("pointer-events","all").on("mouseenter mousemove",function(event,d){ const s=d.summary, fmt=(v)=>transform==="log2"?Math.round(2**v):Math.round(v*100)/100; tooltip(chart.tooltipNode,chart.container,event,d.name,[["Count",String(s.count)],["Q1",String(fmt(s.q1))],["Median",String(fmt(s.median))],["Q3",String(fmt(s.q3))]]); }).on("mouseleave",()=>hideTooltip(chart.tooltipNode));
    }

    function drawBoxes(config) {
      const chart = baseChart(config), plotWidth = WIDTH-MARGIN.left-MARGIN.right, plotHeight=HEIGHT-125-96, transform=config.scaleTransform||"linear";
      titles(chart.svg,config,plotHeight);
      if(config.type==="faceted_boxplot"){
        const facets=list(config.facets), gap=24, panelWidth=(plotWidth-gap*(facets.length-1))/facets.length;
        facets.forEach((facet,i)=>{ const root=chart.svg.append("g").attr("transform",`translate(${MARGIN.left+i*(panelWidth+gap)},125)`); drawBoxPanel(root,facet.groups,panelWidth,plotHeight,transform,chart,i===0); root.append("text").attr("x",panelWidth/2).attr("y",-16).attr("text-anchor","middle").attr("fill",TITLE).style("font-size","13px").style("font-weight",800).text(facet.name); });
      } else {
        const root=chart.svg.append("g").attr("transform",`translate(${MARGIN.left},125)`); drawBoxPanel(root,config.groups,plotWidth,plotHeight,transform,chart,true);
      }
    }

    function drawGeneMatrix(config) {
      const chart=baseChart(config), labels=list(config.sampleLabels), rows=list(config.intersections).map((row)=>({sets:list(row.sets),size:Number(row.size)||0,genes:list(row.genes)})).sort((a,b)=>b.size-a.size||a.sets.length-b.sets.length||a.sets.join("__").localeCompare(b.sets.join("__")));
      if(!labels.length||!rows.length)return;
      const sampleColors=["#d94841","#2563eb","#f59e0b","#0f9d78","#7c3aed"],colorMap=Object.fromEntries(labels.map((label,i)=>[label,sampleColors[i%sampleColors.length]]));
      const sampleSizes=Object.fromEntries(rows.filter((row)=>row.sets.length===1).map((row)=>[row.sets[0],row.size]));
      const blocks=[];for(let i=0;i<rows.length;i+=8)blocks.push(rows.slice(i,i+8));
      const margin={top:52,right:40,bottom:54,left:40},topBarHeight=250,rowHeight=40,columnWidth=132,blockGap=68,blockHeight=topBarHeight+52+labels.length*rowHeight;
      const legendColumns=labels.length>=4?2:Math.min(3,labels.length),legendItemWidth=320,legendTotalWidth=legendColumns*legendItemWidth,maxBlockColumns=d3.max(blocks,(block)=>block.length)||1,maxBlockWidth=Math.max(420,maxBlockColumns*columnWidth);
      const viewportWidth=Math.max(1120,maxBlockWidth+margin.left+76+margin.right,legendTotalWidth+margin.left+margin.right),legendRows=Math.ceil(labels.length/legendColumns),legendTop=margin.top+blocks.length*blockHeight+Math.max(0,blocks.length-1)*blockGap+52,chartHeight=Math.max(700,legendTop+legendRows*68+margin.bottom),globalMax=d3.max(rows,(row)=>row.size)||1;
      chart.svgNode.setAttribute("viewBox",`0 0 ${viewportWidth} ${chartHeight}`);chart.svgNode.setAttribute("preserveAspectRatio","xMidYMid meet");chart.svg.selectAll("*").remove();chart.svg.append("rect").attr("width",viewportWidth).attr("height",chartHeight).attr("fill","#ffffff");
      blocks.forEach((block,bi)=>{const blockTop=margin.top+bi*(blockHeight+blockGap),barBottom=blockTop+topBarHeight,matrixTop=barBottom+52,blockWidth=Math.max(420,block.length*columnWidth),left=(viewportWidth-blockWidth)/2;
        const x=d3.scaleBand().domain(block.map((row)=>row.sets.join("__"))).range([left,left+blockWidth]).paddingInner(.28).paddingOuter(.08),y=d3.scaleLinear().domain([0,globalMax]).nice().range([barBottom,blockTop]);
        chart.svg.append("g").attr("transform",`translate(${left},0)`).call(d3.axisLeft(y).ticks(5).tickFormat((v)=>Number(v).toLocaleString("en-US"))).call((g)=>g.selectAll("text").style("font-size","14px").style("font-weight",700).attr("fill","#4f5b4a")).call((g)=>g.selectAll("path,line").attr("stroke","#74806f").attr("stroke-width",1.2));
        if(bi===0)chart.svg.append("text").attr("x",left+6).attr("y",blockTop-18).attr("fill","#253320").style("font-size","15px").style("font-weight",800).text("Intersection Size");
        chart.svg.append("line").attr("x1",left).attr("x2",left+blockWidth).attr("y1",barBottom).attr("y2",barBottom).attr("stroke","#74806f").attr("stroke-width",1.2);
        block.forEach((row)=>{const id=row.sets.join("__"),x0=x(id),xc=x0+x.bandwidth()/2;
          chart.svg.append("rect").attr("x",x0).attr("y",y(row.size)).attr("width",x.bandwidth()).attr("height",barBottom-y(row.size)).attr("fill","#859b7a").attr("opacity",.88).on("mouseenter mousemove",(event)=>tooltip(chart.tooltipNode,chart.container,event,row.sets.join(" / "),[[row.sets.length===1?"Genes":"Shared Genes",row.size.toLocaleString("en-US")]])).on("mouseleave",()=>hideTooltip(chart.tooltipNode));
          chart.svg.append("text").attr("x",xc).attr("y",y(row.size)-10).attr("text-anchor","middle").attr("fill","#5b4a38").style("font-size",x.bandwidth()>=22?"13px":"11px").style("font-weight",800).text(row.size.toLocaleString("en-US"));
          const active=labels.map((label,i)=>row.sets.includes(label)?matrixTop+i*rowHeight+rowHeight/2:null).filter((v)=>v!==null);if(active.length>1)chart.svg.append("line").attr("x1",xc).attr("x2",xc).attr("y1",d3.min(active)).attr("y2",d3.max(active)).attr("stroke","#7a6655").attr("stroke-width",3).attr("stroke-linecap","round");
          labels.forEach((label,i)=>{const cy=matrixTop+i*rowHeight+rowHeight/2,on=row.sets.includes(label);chart.svg.append("circle").attr("cx",xc).attr("cy",cy).attr("r",7.5).attr("fill",on?colorMap[label]:"#d9ded3").attr("stroke",on?"#ffffff":"#b8c2b1").attr("stroke-width",on?2.2:1.1);});
        });
      });
      const legendLeft=(viewportWidth-legendTotalWidth)/2+72;labels.forEach((label,i)=>{const col=i%legendColumns,row=Math.floor(i/legendColumns),x=legendLeft+col*legendItemWidth,y=legendTop+row*68,color=colorMap[label];chart.svg.append("rect").attr("x",x).attr("y",y+2).attr("width",58).attr("height",16).attr("rx",5).attr("fill",color).attr("opacity",.9);chart.svg.append("text").attr("x",x+70).attr("y",y+4).attr("dominant-baseline","hanging").attr("fill","#253320").style("font-size","13px").style("font-weight",700).text(label);chart.svg.append("text").attr("x",x+70).attr("y",y+30).attr("dominant-baseline","hanging").attr("fill","#5b4a38").style("font-size","13px").style("font-weight",800).text((sampleSizes[label]||0).toLocaleString("en-US"));});
    }    function drawSiteProfile(config) {
      const chart = baseChart(config), panels = list(config.panels);
      if (!panels.length) return;
      const panel = panels[0], width = WIDTH, height = HEIGHT, axis = AXIS;
      chart.svgNode.setAttribute("viewBox", `0 0 ${width} ${height}`);
      chart.svg.selectAll("*").remove();
      chart.svg.append("rect").attr("width", width).attr("height", height).attr("fill", "#fffcf6");
      chart.svg.append("text").attr("x",112).attr("y",34).attr("fill",TITLE).style("font-size","15px").style("font-weight",800).text(panel.title||"");
      if (panel.type === "density") {
        const margin={top:64,right:128,bottom:130,left:112}, plotWidth=width-margin.left-margin.right, plotHeight=height-margin.top-margin.bottom;
        const series=list(panel.series), points=series.flatMap((s)=>list(s.values));
        const x=d3.scaleLinear().domain(panel.xDomain||d3.extent(points,(d)=>Number(d.x))).range([0,plotWidth]);
        const yDomain=panel.yDomain||[0,Math.max(1e-12,d3.max(points,(d)=>Number(d.density))||0)*1.08];
        const y=d3.scaleLinear().domain(yDomain).nice().range([plotHeight,0]);
        const root=chart.svg.append("g").attr("transform",`translate(${margin.left},${margin.top})`), yTicks=list(panel.yTicks).length?list(panel.yTicks):y.ticks(5);
        root.selectAll(".grid").data(yTicks).join("line").attr("x2",plotWidth).attr("y1",y).attr("y2",y).attr("stroke",GRID).attr("stroke-width",1);
        root.append("g").call(d3.axisLeft(y).tickValues(yTicks).tickSize(6).tickPadding(10)).call((g)=>g.select(".domain").remove()).call((g)=>g.selectAll(".tick line").attr("stroke",axis).attr("stroke-width",1.2)).call((g)=>g.selectAll("text").attr("fill",axis).style("font-size","12px").style("font-weight",700));
        root.append("line").attr("y2",plotHeight).attr("stroke",axis).attr("stroke-width",1.2);
        root.append("g").attr("transform",`translate(0,${plotHeight})`).call(d3.axisBottom(x).ticks(6).tickSize(0).tickPadding(12)).call((g)=>g.select(".domain").attr("stroke",axis).attr("stroke-width",1.2)).call((g)=>g.selectAll("text").attr("fill",axis).style("font-size","12px").style("font-weight",700));
        root.append("text").attr("x",plotWidth/2).attr("y",plotHeight+52).attr("text-anchor","middle").attr("fill",axis).style("font-size","15px").style("font-weight",700).text(panel.xLabel||"");
        root.append("text").attr("transform",`translate(-78,${plotHeight/2}) rotate(-90)`).attr("text-anchor","middle").attr("fill",axis).style("font-size","15px").style("font-weight",700).text(panel.yLabel||"");
        list(panel.guideLines).forEach((value)=>root.append("line").attr("x1",x(Number(value))).attr("x2",x(Number(value))).attr("y2",plotHeight).attr("stroke",axis).attr("stroke-width",1.2).attr("stroke-dasharray","3,6"));
        const line=d3.line().x((d)=>x(Number(d.x))).y((d)=>y(Number(d.density))).curve(d3.curveMonotoneX);
        series.forEach((entry)=>root.append("path").datum(list(entry.values)).attr("fill","none").attr("stroke",entry.color).attr("stroke-width",3.6).attr("d",line).on("mouseenter mousemove",(event)=>{const target=x.invert(d3.pointer(event,root.node())[0]);const values=list(entry.values);const datum=values.reduce((best,value)=>Math.abs(Number(value.x)-target)<Math.abs(Number(best.x)-target)?value:best,values[0]);tooltip(chart.tooltipNode,chart.container,event,entry.name,[["Panel",panel.title||""],["Position",Number(datum.x).toFixed(1)],["Density",Number(datum.density).toFixed(4)]]);}).on("mouseleave",()=>hideTooltip(chart.tooltipNode)));
        if(series.length>1) legend(chart.svg,series,height-30);
      } else if (panel.type === "heatmap") {
        const margin={top:64,right:128,bottom:74,left:112}, plotWidth=width-margin.left-margin.right, plotHeight=height-margin.top-margin.bottom;
        const rows=list(panel.values).map(list), nr=rows.length, nc=nr?rows[0].length:0, values=rows.flat().map(Number).sort(d3.ascending);
        const maxValue=d3.quantile(values,Number(panel.colorMaxQuantile)||.92)||d3.max(values)||1;
        const color=d3.scaleSequential().domain([0,maxValue]).interpolator(d3.interpolateRgbBasis(["#ffffff","#f7cf7a","#ee9b4b","#da6a33","#b23a2c","#7f1d1d"]));
        const root=chart.svg.append("g").attr("transform",`translate(${margin.left},${margin.top})`), cw=plotWidth/Math.max(1,nc), rh=plotHeight/Math.max(1,nr), cells=[];
        rows.forEach((row,r)=>row.forEach((value,c)=>cells.push({r,c,value:Number(value)})));
        root.append("rect").attr("width",plotWidth).attr("height",plotHeight).attr("fill",panel.backgroundColor||"#ffffff");
        root.selectAll("rect.heat-cell").data(cells).join("rect").attr("class","heat-cell").attr("x",(d)=>d.c*cw).attr("y",(d)=>d.r*rh).attr("width",cw+.2).attr("height",rh+.2).attr("fill",(d)=>color(Math.min(d.value,maxValue)));
        const y=d3.scaleLinear().domain([1,Math.max(1,nr)]).range([rh/2,plotHeight-rh/2]);
        root.append("g").call(d3.axisLeft(y).tickValues(Array.from(new Set([1,Math.max(1,Math.round(nr/2)),Math.max(1,nr)]))).tickSize(8).tickPadding(8)).call((g)=>g.select(".domain").remove()).call((g)=>g.selectAll(".tick line").attr("stroke",axis).attr("stroke-width",1.2)).call((g)=>g.selectAll("text").attr("fill",axis).style("font-size","12px").style("font-weight",700));
        root.append("line").attr("y2",plotHeight).attr("stroke",axis).attr("stroke-width",1.2);
        const xd=panel.xDomain||[panel.xMin,panel.xMax], x=d3.scaleLinear().domain(xd).range([0,plotWidth]);
        root.append("g").attr("transform",`translate(0,${plotHeight})`).call(d3.axisBottom(x).ticks(5).tickSize(8).tickPadding(8)).call((g)=>g.select(".domain").attr("stroke",axis).attr("stroke-width",1.2)).call((g)=>g.selectAll(".tick line").attr("stroke",axis).attr("stroke-width",1.2)).call((g)=>g.selectAll("text").attr("fill",axis).style("font-size","12px").style("font-weight",700));
      }
    }    function drawMeta(config) {
      const chart = baseChart(config), series = list(config.series);
      if (!series.length) return;
      const metaHeight = 520, margin = { top: 58, right: 28, bottom: 112, left: 76 };
      const plotWidth = WIDTH - margin.left - margin.right, plotHeight = metaHeight - margin.top - margin.bottom;
      chart.svgNode.setAttribute("viewBox", `0 0 ${WIDTH} ${metaHeight}`);
      const svg = chart.svg, panel = svg.append("g");
      const x = d3.scaleLinear().domain(config.xDomain || [0, 1]).range([margin.left, WIDTH - margin.right]);
      const y = d3.scaleLinear().domain(config.yDomain || [0, 1]).range([margin.top + plotHeight, margin.top]);
      panel.append("rect").attr("x", margin.left).attr("y", margin.top).attr("width", plotWidth).attr("height", plotHeight).attr("fill", "#EBEBEB");
      panel.selectAll(".meta-grid-x").data(d3.range(0, 1.0001, .25)).join("line").attr("x1", x).attr("x2", x).attr("y1", margin.top).attr("y2", margin.top + plotHeight).attr("stroke", "rgba(255,255,255,.95)").attr("stroke-width", 1.1);
      const yTicks = y.ticks(5);
      panel.selectAll(".meta-grid-y").data(yTicks).join("line").attr("x1", margin.left).attr("x2", WIDTH-margin.right).attr("y1", y).attr("y2", y).attr("stroke", "rgba(255,255,255,.95)").attr("stroke-width", 1.1);
      panel.append("g").attr("transform", `translate(${margin.left},0)`).call(d3.axisLeft(y).tickValues(yTicks).tickSize(6).tickPadding(8)).call((g)=>g.select(".domain").remove()).call((g)=>g.selectAll("text").attr("fill","#4c4c4c").style("font-size","14px"));
      panel.append("line").attr("x1",margin.left).attr("x2",margin.left).attr("y1",margin.top).attr("y2",margin.top+plotHeight).attr("stroke","#4c4c4c").attr("stroke-width",1.5);
      svg.append("text").attr("x",margin.left).attr("y",50).attr("fill","#111").style("font-size","15px").style("font-weight",700).text(config.title||"");
      svg.append("text").attr("transform","translate(30,250) rotate(-90)").attr("fill","#111").style("font-size","15px").style("font-weight",700).text(config.yLabel||"");
      const area=d3.area().x((d)=>x(Number(d.x))).y0(y(0)).y1((d)=>y(Number(d.density))).curve(d3.curveBasis);
      const line=d3.line().x((d)=>x(Number(d.x))).y((d)=>y(Number(d.density))).curve(d3.curveBasis);
      const ciArea=d3.area().x((d)=>x(Number(d.x))).y0((d)=>y(Number(d.confidenceDown))).y1((d)=>y(Number(d.confidenceUp))).curve(d3.curveBasis);
      series.forEach((entry)=>{
        const values=list(entry.values);
        panel.append("path").datum(values).attr("fill",entry.color).attr("opacity",.18).attr("d",area);
        if(config.showCI && values.every((d)=>Number.isFinite(Number(d.confidenceDown))&&Number.isFinite(Number(d.confidenceUp)))) panel.append("path").datum(values).attr("fill",entry.color).attr("opacity",.12).attr("d",ciArea);
        panel.append("path").datum(values).attr("fill","none").attr("stroke",entry.color).attr("stroke-width",4).attr("d",line);
        panel.append("path").datum(values).attr("fill","none").attr("stroke","transparent").attr("stroke-width",18).attr("d",line).style("pointer-events","stroke")
          .on("mouseenter mousemove",(event)=>{const target=x.invert(d3.pointer(event,panel.node())[0]);const datum=values.reduce((best,value)=>Math.abs(Number(value.x)-target)<Math.abs(Number(best.x)-target)?value:best,values[0]);tooltip(chart.tooltipNode,chart.container,event,entry.name,[["Position",Number(datum.x).toFixed(3)],[config.yLabel||"Frequency",Number(datum.density).toFixed(3)]]);}).on("mouseleave",()=>hideTooltip(chart.tooltipNode));
      });
      const components=config.components||{}, componentBase=Number(components.trackBase??-.03), labelY=Number(components.labelY??-.08);
      list(components.separators).forEach((separator)=>panel.append("line").attr("x1",x(Number(separator.x))).attr("x2",x(Number(separator.x))).attr("y1",y(Number(separator.y1))).attr("y2",y(Number(separator.y2))).attr("stroke","#111").attr("stroke-width",1.6).attr("stroke-dasharray","2,6"));
      const minimum={promoter:2,tail:2,utr5:4,utr3:4,cds:6,ncrna:6,rna:6};
      list(components.segments).forEach((segment)=>{const center=y(componentBase),height=Math.max(minimum[segment.component]||3,Math.abs(y(componentBase-Number(segment.height))-y(componentBase+Number(segment.height)))),gx=x(Number(segment.start)),gw=Math.max(1,x(Number(segment.end))-gx);panel.append("rect").attr("x",gx).attr("y",Math.round(center-height/2)+.5).attr("width",gw).attr("height",Math.max(1,Math.round(height))).attr("fill","#fff").attr("fill-opacity",Number(segment.alpha??1)).attr("stroke","#111").attr("stroke-width",1.1);if(segment.label)panel.append("text").attr("x",x(Number(segment.mid))).attr("y",y(labelY)+12).attr("text-anchor","middle").attr("fill","#222").style("font-size","14px").style("font-weight",700).text(segment.label);});
      legend(svg,series,metaHeight-85);
    }
    function draw(config) {
      if (["grouped_bar","stacked_bar"].includes(config.type)) drawBars(config);
      else if (["boxplot","faceted_boxplot"].includes(config.type)) drawBoxes(config);
      else if (config.type === "gene_matrix") drawGeneMatrix(config);
      else if (config.type === "site_profile") drawSiteProfile(config);
      else if (config.type === "meta_plot") drawMeta(config);
      else baseChart(config);
    }
    return { renderValue: (x) => draw(x.config || {}), resize: function() {} };
  }
});
